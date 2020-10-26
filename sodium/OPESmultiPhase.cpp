/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2020 of Michele Invernizzi.

The opes module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The opes module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "bias/Bias.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "core/Atoms.h"
#include "tools/Communicator.h"
#include "tools/File.h"

namespace PLMD {
namespace opes {

//+PLUMEDOC BIAS OPES_MULTIPHASE
/*
\par Examples

OPES_MULTIPHASE ...
  LABEL=test
  ARG=ene,vol
  PACE=50
  TEMP=300
  MIN_TEMP=260
  MAX_TEMP=350
  PRESSURE=0.06022140857
  MAX_PRESSURE=180.66422571
  MIN_PRESSURE=0.06022140857
  SIGMA=0.18
  MIN_CV=-1
  MAX_CV=1
... OPES_MULTIPHASE


*/
//+ENDPLUMEDOC

class OPESmultiPhase : public bias::Bias {

private:
  bool isFirstStep_;
  bool afterCalculate_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  unsigned long counter_;

  unsigned stride_;
  unsigned obs_steps_;
  std::vector<double> obs_ene_;
  std::vector<double> obs_vol_;
  std::vector<double> obs_cv_;
  unsigned steps_beta_;
  unsigned steps_pres_;
  unsigned tot_umbrellas_;
  unsigned P0_contribution_;

  double beta_;
  double pres_;
  double sigma_;
  std::vector<double> beta_p_;
  std::vector<double> pres_p_;
  std::vector<double> center_;
  std::vector< std::vector< std::vector<double> > > deltaF_;
  double rct_;
  double my_rct_;
  double barrier_;
  double current_bias_;

  double border_weight_;
  double tot_steps_;

  bool calc_work_;
  double work_;
  std::vector< std::vector< std::vector<double> > > old_deltaF_;

  OFile deltaFsOfile_;
  unsigned print_stride_;

public:
  OPESmultiPhase(const ActionOptions&);
  inline long double get_weight(const double,const double,const double,const double,const double,const double,const double) const;
  void calculate() override;
  void update() override;
  void init_integration_grid();
  void init_from_obs();
  unsigned estimate_steps(const double,const double,const std::vector<double>&,const std::string) const;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(OPESmultiPhase,"OPES_MULTIPHASE")

void OPESmultiPhase::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","PACE","10","how often the bias is updated");
  keys.add("compulsory","OBSERVATION_STEPS","100","number of unbiased initial steps to collect statistics."
           " When using STEPS_TEMP and STEPS_PRESSURE it is forced to 1");
//tempertature stuff
  keys.add("compulsory","TEMP","-1","temperature. If not specified tries to get it from MD engine");
  keys.add("compulsory","MIN_TEMP","the minimum of the temperature range");
  keys.add("compulsory","MAX_TEMP","the maximum of the temperature range");
  keys.add("optional","STEPS_TEMP","manually set the number of intermediate temperatures");
//pressure stuff
  keys.add("compulsory","PRESSURE","simulation pressure");
  keys.add("compulsory","MIN_PRESSURE","the minimum of the pressure range");
  keys.add("compulsory","MAX_PRESSURE","the maximum of the pressure range");
  keys.add("optional","STEPS_PRESSURE","manually set the number of intermediate pressures");
//umbrella stuff
  keys.add("compulsory","SIGMA","sigma of the umbrella Gaussians");
  keys.add("compulsory","MIN_CV","the minimum of the CV range to be explored");
  keys.add("compulsory","MAX_CV","the maximum of the CV range to be explored");
  keys.add("optional","BARRIER_CV","a guess of the free energy barrier to be overcome (better to stay higher than lower)");
  keys.addFlag("ADD_P0",false,"add also P0 to the target distribution, useful to make sure the target is broader than P0");
//deltaFs file
  keys.add("compulsory","FILE","DELTAFS","a file with the estimate of the relative \\f$\\Delta F\\f$ for each component of the target");
  keys.add("optional","PRINT_STRIDE","( default=100 ) stride for printing to DELTAFS file");
  keys.add("optional","FMT","specify format for DELTAFS file");
//miscellaneous
  keys.add("optional","BORDER_WEIGHT","set it greater than 1 to obtain better sampling of the max and min thermodynamics conditions");
  keys.addFlag("CALC_WORK",false,"calculate the work done by the bias between each update");
  keys.addFlag("WALKERS_MPI",false,"switch on MPI version of multiple walkers");
  keys.addFlag("SERIAL",false,"perform calculations in serial");
  keys.use("RESTART");

//output components
  componentsAreNotOptional(keys);
  keys.addOutputComponent("rct","WALKERS_MPI","single walker estimate of \\f$c(t)\\f$: \\f$\\frac{1}{\\beta}\\log \\lange e^{\\beta V} \\rangle\\f$, all walkers should converge to the same value");
  keys.addOutputComponent("work","CALC_WORK","work done by the bias between each update"); //calculating this maybe is only a useless overhead...
}

OPESmultiPhase::OPESmultiPhase(const ActionOptions&ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , counter_(0)
  , rct_(0)
  , my_rct_(0)
  , work_(0)
  , print_stride_(100)
{
  //TODO is there a way to check that ARG is actually energy,volume,CV?
  plumed_massert(getNumberOfArguments()==3,"only ENERGY, VOLUME and a CV should be given as ARG");

//set beta_
  const double Kb=plumed.getAtoms().getKBoltzmann();
  double KbT=plumed.getAtoms().getKbT();
  double temp=-1;
  parse("TEMP",temp);
  if(temp>0)
  {
    if(KbT>0 && std::abs(KbT-Kb*temp)>1e-4)
      log.printf(" +++ WARNING +++ using TEMP=%g while MD engine uses %g\n",temp,KbT/Kb);
    KbT=Kb*temp;
  }
  plumed_massert(KbT>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
  beta_=1./KbT;
  temp=1./(Kb*beta_);

//set temp range
  double min_temp;
  double max_temp;
  steps_beta_=0;
  parse("MIN_TEMP",min_temp);
  parse("MAX_TEMP",max_temp);
  parse("STEPS_TEMP",steps_beta_);
  plumed_massert(max_temp>=min_temp,"MAX_TEMP cannot be smaller than MIN_TEMP");
  const double tol=1e-3; //if temp is taken from MD engine it might be numerically slightly different
  if(temp<(1-tol)*min_temp || temp>(1+tol)*max_temp)
    log.printf(" +++ WARNING +++ running at TEMP=%g which is outside the chosen temperature range\n",temp);
  if(min_temp==max_temp)
  {
    plumed_massert(std::abs((temp-min_temp)/min_temp)<tol,"if MIN_TEMP = MAX_TEMP, they should be equal to TEMP");
    plumed_massert(steps_beta_==0 || steps_beta_==1,"cannot have multiple steps when MIN_TEMP = MAX_TEMP");
    steps_beta_=1;
  }
  beta_p_.resize(2);
  beta_p_[0]=1./(Kb*max_temp); //min_beta
  beta_p_[1]=1./(Kb*min_temp); //max_beta

//set pres
  steps_pres_=0;
  parse("PRESSURE",pres_);
  pres_p_.resize(2);
  parse("MIN_PRESSURE",pres_p_[0]);
  parse("MAX_PRESSURE",pres_p_[1]);
  parse("STEPS_PRESSURE",steps_pres_);
  plumed_massert(pres_p_[1]>=pres_p_[0],"MAX_PRESSURE is smaller than MIN_PRESSURE");
  if(pres_p_[0]==pres_p_[1])
  {
    plumed_massert(pres_==pres_p_[0],"if MIN_PRESSURE = MAX_PRESSURE, they should be equal to PRESSURE");
    plumed_massert(steps_pres_==0 || steps_pres_==1,"cannot have multiple steps when MIN_PRESSURE = MAX_PRESSURE");
    steps_pres_=1;
  }

//set other stuff
  parse("PACE",stride_);
  parse("OBSERVATION_STEPS",obs_steps_);
  plumed_massert(obs_steps_!=0,"minimum is OBSERVATION_STEPS=1");
  if(steps_beta_!=0 && steps_pres_!=0)
    obs_steps_=1;
  obs_ene_.resize(obs_steps_);
  obs_vol_.resize(obs_steps_);
  obs_cv_.resize(obs_steps_);

  border_weight_=1;
  parse("BORDER_WEIGHT",border_weight_);

//set umbrellas
  parse("SIGMA",sigma_);
  double min_cv;
  double max_cv;
  parse("MIN_CV",min_cv);
  parse("MAX_CV",max_cv);
  plumed_massert(min_cv<max_cv,"MIN_CV should be smaller than MAX_CV");
  tot_umbrellas_=1+std::round((max_cv-min_cv)/sigma_);
  center_.resize(tot_umbrellas_);
  for(unsigned k=0; k<tot_umbrellas_; k++)
    center_[k]=(max_cv-min_cv)/(tot_umbrellas_-1)*k+min_cv;
  barrier_=std::numeric_limits<double>::infinity();
  parse("BARRIER_CV",barrier_);
  bool add_P0=false;
  parseFlag("ADD_P0",add_P0);
  if(add_P0)
    P0_contribution_=1;
  else
    P0_contribution_=0;

//deltaFs file
  std::string deltaFsFileName;
  parse("FILE",deltaFsFileName);
  parse("PRINT_STRIDE",print_stride_);
  std::string fmt;
  parse("FMT",fmt);

//work flag
  parseFlag("CALC_WORK",calc_work_);

//multiple walkers //TODO implement also external mw for cp2k
  unsigned walker_rank;
  bool walkers_mpi=false;
  parseFlag("WALKERS_MPI",walkers_mpi);
  if(walkers_mpi)
  {
    if(comm.Get_rank()==0)//multi_sim_comm works on first rank only
    {
      NumWalkers_=multi_sim_comm.Get_size();
      walker_rank=multi_sim_comm.Get_rank();
    }
    comm.Bcast(NumWalkers_,0); //if each walker has more than one processor update them all
    comm.Bcast(walker_rank,0);
  }
  else
  {
    NumWalkers_=1;
    walker_rank=0;
  }

//parallelization stuff
  NumParallel_=comm.Get_size();
  rank_=comm.Get_rank();
  bool serial=false;
  parseFlag("SERIAL",serial);
  if(serial)
  {
    log.printf(" -- SERIAL: running without loop parallelization\n");
    NumParallel_=1;
    rank_=0;
  }

  checkRead();

//printing some info
  log.printf("  Beta = %g\n",beta_);
  log.printf("  Running at TEMP = %g\n",1./(Kb*beta_));
  log.printf("  Running at PRESSURE = %g\n",pres_);
  log.printf("  Updating the bias with PACE = %u\n",stride_);
  log.printf("  Initial unbiased observation done for OBSERVATION_STEPS = %u\n",obs_steps_);
  if(steps_beta_!=0 && steps_pres_!=0)
    log.printf(" +++ WARNING +++ STEPS_TEMP and STEPS_PRES are used, thus OBSERVATION_STEP is set to 1\n"
               "                 Reference energy and initial deltaFs can be set with a fake RESTART\n");
  if(barrier_!=std::numeric_limits<double>::infinity())
    log.printf("  guess for free energy BARRIER_CV = %g\n",barrier_);
  if(add_P0)
    log.printf(" --- ADD_P0: adding P0 to the target\n");
  if(walkers_mpi)
    log.printf("  WALKERS_MPI: multiple walkers will share the same bias via mpi\n");
  if(NumWalkers_>1)
  {
    log.printf("  Using multiple walkers\n");
    log.printf("    number of walkers: %u\n",NumWalkers_);
    log.printf("    walker rank: %u\n",walker_rank);
  }
  if(NumParallel_>1)
    log.printf("  Using multiple threads per simulation: %u\n",NumParallel_);

//restart if needed
  if(getRestart())
  {
    IFile ifile;
    ifile.link(*this);
    if(ifile.FileExist(deltaFsFileName))
    {
      log.printf("  RESTART from: %s\n",deltaFsFileName.c_str());
      log.printf(" +++ make sure all simulation options are consistent! +++\n");
      ifile.open(deltaFsFileName);
      std::vector<std::string> deltaFsName;
      ifile.scanFieldList(deltaFsName);//TODO add a check to see if the temp/pres grid is the same?
      const unsigned pre_f=2; //time and rct
      const unsigned post_f=1; //print_stride
      plumed_massert(deltaFsName.size()-pre_f-post_f>=2,"RESTART - fewer than expected FIELDS found in '"+deltaFsFileName+"' file");
//      const std::string lastW=deltaFsName[deltaFsName.size()-post_f-1];
//      const std::size_t first=lastW.find_first_of("_");
//      const std::size_t last=lastW.find_last_of("_");
//      plumed_massert(first!=last,"RESTART - something wrong with '"+deltaFsFileName+"' file: could not find 2 steps in "+lastW);
      std::string lastW=deltaFsName[deltaFsName.size()-post_f-1];
      try
      {
        std::size_t last=lastW.find_last_of("_");
        unsigned used_umbrellas=std::stoi(lastW.substr(last+1,lastW.size()-last-1));
        plumed_massert(used_umbrellas==tot_umbrellas_,"mismatch in the number of umbrella steps");
        lastW=lastW.substr(0,last);
        last=lastW.find_last_of("_");
        steps_pres_=std::stoi(lastW.substr(last+1,lastW.size()-last-1));
        lastW=lastW.substr(0,last);
        last=lastW.find_last_of("_");
        steps_beta_=std::stoi(lastW.substr(last+1,lastW.size()-last-1));
      }
      catch(std::exception const & e)
      {
        error("RESTART - something wrong with '"+deltaFsFileName+"' file, expected FIELDS format is: deltaF_b_p_u");
      }
      plumed_massert(steps_beta_*steps_pres_*tot_umbrellas_==deltaFsName.size()-pre_f-post_f,"RESTART - something wrong with '"+deltaFsFileName+"' file: steps not matching size");
      init_integration_grid();
      deltaF_.resize(steps_beta_,std::vector< std::vector<double> >(steps_pres_,std::vector<double>(tot_umbrellas_)));
      obs_steps_=0; //avoid initializing again
      //read steps from file
      int restart_stride;
      ifile.scanField("print_stride",restart_stride);
      plumed_massert(restart_stride==(int)print_stride_,"also PRINT_STRIDE must be consistent to avoid problems with multiple restarts");
      ifile.allowIgnoredFields(); //this allows for multiple restart, but without checking for consistency between them!
      double time;
      while(ifile.scanField("time",time)) //room for improvements: only last line is important
      {
        if(calc_work_)
          old_deltaF_=deltaF_;
        ifile.scanField("rct",rct_);
        for(unsigned i=0; i<steps_beta_; i++)
          for(unsigned j=0; j<steps_pres_; j++)
            for(unsigned k=0; k<tot_umbrellas_; k++)
              ifile.scanField(deltaFsName[pre_f+i*steps_pres_*tot_umbrellas_+j*tot_umbrellas_+k],deltaF_[i][j][k]);
        ifile.scanField();
        counter_++;
      }
      log.printf("  Successfully read %d lines, up to t=%g\n",counter_,time);
      counter_=(1+(counter_-1)*print_stride_)*NumWalkers_; //adjust counter
      if(NumWalkers_>1)
      {
        my_rct_=rct_; //better than 0
        log.printf(" +++ WARNING +++ the single walker rct estimate is resetted to the global rct\n");
      }
      ifile.reset(false);
      ifile.close();
    }
    else
      log.printf(" +++ WARNING +++ restart requested, but no '%s' file found!\n",deltaFsFileName.c_str());
  }
//sync all walkers to avoid opening files before reding is over (see also METAD)
  comm.Barrier();
  if(comm.Get_rank()==0 && walkers_mpi)
    multi_sim_comm.Barrier();

//setup deltaFs file, without opening it
  deltaFsOfile_.link(*this);
  if(NumWalkers_>1)
  {
    if(walker_rank>0)
      deltaFsFileName="/dev/null"; //only first walker writes on file
    deltaFsOfile_.enforceSuffix("");
  }
  deltaFsOfile_.open(deltaFsFileName);
  if(fmt.length()>0)
    deltaFsOfile_.fmtField(" "+fmt);
  deltaFsOfile_.setHeavyFlush(); //do I need it?
  deltaFsOfile_.addConstantField("print_stride");
  deltaFsOfile_.printField("print_stride",(int)print_stride_);

//add output components
  if(NumWalkers_>1)
  {
    addComponent("rct");
    componentIsNotPeriodic("rct");
    getPntrToComponent("rct")->set(my_rct_);
  }
  if(calc_work_)
  {
    addComponent("work");
    componentIsNotPeriodic("work");
  }

//Bibliography
  log.printf("  Bibliography: ");
  log<<plumed.cite("M. Invernizzi, P.M. Piaggi, and M. Parrinello, arXiv:2007.03055 (2020)");
  log.printf("\n");
}

inline long double OPESmultiPhase::get_weight(const double beta_p,const double ene,const double pres_p,const double vol,const double center, const double cv, const double shift) const
{ //FIXME avoid long double?
  return std::exp(static_cast<long double>((beta_-beta_p)*ene+(beta_*pres_-beta_p*pres_p)*vol-0.5*std::pow((center-cv)/sigma_,2)+beta_*shift));
}

void OPESmultiPhase::calculate()
{
  if(obs_steps_>0) //no bias before initialization
    return;

  const double ene=getArgument(0);
  const double vol=getArgument(1);
  const double cv=getArgument(2);
  long double sum=0;
  long double der_sum_ene=0;
  long double der_sum_vol=0;
  long double der_sum_cv=0;
  for(unsigned i=rank_; i<steps_beta_; i+=NumParallel_)
  {
    for(unsigned j=0; j<steps_pres_; j++)
    {
      double w=1;
      if(i==0 || i==steps_beta_-1 || j==0 || j==steps_pres_-1)
        w=border_weight_;
      for(unsigned k=0; k<tot_umbrellas_; k++)
      {
        const long double add_ijk=w*get_weight(beta_p_[i],ene,pres_p_[j],vol,center_[k],cv,deltaF_[i][j][k]);
        sum+=add_ijk;
        der_sum_ene+=(beta_-beta_p_[i])*add_ijk;
        der_sum_vol+=(beta_*pres_-beta_p_[i]*pres_p_[j])*add_ijk;
        der_sum_cv-=(cv-center_[k])/sigma_/sigma_*add_ijk;
      }
    }
  }
  if(NumParallel_>1)
  {
    comm.Sum(sum);
    comm.Sum(der_sum_ene);
    comm.Sum(der_sum_vol);
    comm.Sum(der_sum_cv);
  }
  sum+=P0_contribution_;

  current_bias_=-1./beta_*std::log(sum/tot_steps_);
  setBias(current_bias_);

  setOutputForce(0,1./beta_*der_sum_ene/sum);
  setOutputForce(1,1./beta_*der_sum_vol/sum);
  setOutputForce(2,1./beta_*der_sum_cv/sum);

//calculate work
  if(calc_work_)
  {
    long double old_sum=0;
    for(unsigned i=rank_; i<steps_beta_; i+=NumParallel_)
    {
      for(unsigned j=0; j<steps_pres_; j++)
      {
        double w=1;
        if(i==0 || i==steps_beta_-1 || j==0 || j==steps_pres_-1)
          w=border_weight_;
        for(unsigned k=0; k<tot_umbrellas_; k++)
          old_sum+=w*get_weight(beta_p_[i],ene,pres_p_[j],vol,center_[k],cv,old_deltaF_[i][j][k]);
      }
    }
    if(NumParallel_>1)
      comm.Sum(old_sum);
    old_sum+=P0_contribution_;
    work_+=-1./beta_*std::log(sum/old_sum);
  }

  afterCalculate_=true;
}

void OPESmultiPhase::update()
{
  if(getStep()%stride_!=0)
    return;
  if(isFirstStep_) //skip very first step, as in METAD
  {
    isFirstStep_=false;
    if(obs_steps_!=1) //if obs_steps_==1 go on with initialization
      return;
  }
  if(obs_steps_>0)
  {
    obs_ene_[counter_]=getArgument(0);
    obs_vol_[counter_]=getArgument(1);
    obs_cv_[counter_]=getArgument(2);
    counter_++;
    if(counter_==obs_steps_)
    {
      log.printf("\nAction OPES_MULTIPHASE\n");
      init_from_obs();
      log.printf("Finished initialization\n\n");
      counter_=NumWalkers_; //all preliminary observations count 1
      obs_steps_=0; //no more observation
    }
    return;
  }
  plumed_massert(afterCalculate_,"OPESmultiPhase::update() must be called after OPESmultiPhase::calculate() to work properly");
  afterCalculate_=false;

//work done by the bias in one iteration
  if(calc_work_)
  {
    getPntrToComponent("work")->set(work_);
    work_=0;
    old_deltaF_=deltaF_;
  }

//update averages. This assumes that calculate() always runs before update(), thus uses current_bias_
  const double ene=getArgument(0);
  const double vol=getArgument(1);
  const double cv=getArgument(2);
  if(NumWalkers_==1)
  { //log1p(arg)=log(1+arg)
    counter_++;
    const double increment=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(current_bias_-rct_)))/(counter_-1.));
    for(unsigned i=0; i<steps_beta_; i++)//TODO is more likely that tot_umbrellas>steps_beta and steps_pres, so why not put it first and parallelize it?
      for(unsigned j=0; j<steps_pres_; j++)
        for(unsigned k=0; k<tot_umbrellas_; k++)
          deltaF_[i][j][k]+=increment-1./beta_*std::log1p(get_weight(beta_p_[i],ene,pres_p_[j],vol,center_[k],cv,current_bias_-rct_+deltaF_[i][j][k])/(counter_-1.));
    rct_+=increment+1./beta_*std::log1p(-1./counter_);
  }
  else
  {
    std::vector<double> all_bias(NumWalkers_);
    std::vector<double> all_ene(NumWalkers_);
    std::vector<double> all_vol(NumWalkers_);
    std::vector<double> all_cv(NumWalkers_);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(current_bias_,all_bias);
      multi_sim_comm.Allgather(ene,all_ene);
      multi_sim_comm.Allgather(vol,all_vol);
      multi_sim_comm.Allgather(cv,all_cv);
    }
    comm.Bcast(all_bias,0);
    comm.Bcast(all_ene,0);
    comm.Bcast(all_vol,0);
    comm.Bcast(all_cv,0);
    for(unsigned w=0; w<NumWalkers_; w++)
    {
      counter_++;
      const double increment_w=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(all_bias[w]-rct_)))/(counter_-1.));
      for(unsigned i=0; i<steps_beta_; i++)
        for(unsigned j=0; j<steps_pres_; j++)
          for(unsigned k=0; k<tot_umbrellas_; k++)
            deltaF_[i][j][k]+=increment_w-1./beta_*std::log1p(get_weight(beta_p_[i],all_ene[w],pres_p_[j],all_vol[w],center_[k],all_cv[w],all_bias[w]-rct_+deltaF_[i][j][k])/(counter_-1.));
      rct_+=increment_w+1./beta_*std::log1p(-1./counter_);
    }
    //calc single walker rct
    const unsigned single_counter=counter_/NumWalkers_;
    const double increment=1./beta_*std::log1p(std::exp(static_cast<long double>(beta_*(current_bias_-my_rct_)))/(single_counter-1.));
    my_rct_+=increment+1./beta_*std::log1p(-1./single_counter);
    getPntrToComponent("rct")->set(my_rct_);
  }

//write to file
  if((counter_/NumWalkers_-1)%print_stride_==0)
  {
    deltaFsOfile_.printField("time",getTime());
    deltaFsOfile_.printField("rct",rct_);
    for(unsigned i=0; i<steps_beta_; i++)
      for(unsigned j=0; j<steps_pres_; j++)
        for(unsigned k=0; k<tot_umbrellas_; k++)
          deltaFsOfile_.printField("deltaF_"+std::to_string(i+1)+"_"+std::to_string(j+1)+"_"+std::to_string(k+1),deltaF_[i][j][k]);
    deltaFsOfile_.printField();
  }
}

void OPESmultiPhase::init_integration_grid()
{
  plumed_massert(steps_beta_*steps_pres_*tot_umbrellas_>0,"zero steps given for integration grid");
//initialize beta grid
  const double min_beta=beta_p_[0];
  const double max_beta=beta_p_[1];
  beta_p_.clear();
  beta_p_.resize(steps_beta_);
  if(steps_beta_==1)
    beta_p_[0]=(min_beta+max_beta)/2.;
  else
    for(unsigned i=0; i<steps_beta_; i++)
      beta_p_[i]=min_beta+i*(max_beta-min_beta)/(steps_beta_-1);
//initialize pres grid
  const double min_pres=pres_p_[0];
  const double max_pres=pres_p_[1];
  pres_p_.clear();
  pres_p_.resize(steps_pres_);
  if(steps_pres_==1)
    pres_p_[0]=(min_pres+max_pres)/2.;
  else
    for(unsigned j=0; j<steps_pres_; j++)
      pres_p_[j]=min_pres+j*(max_pres-min_pres)/(steps_pres_-1);

//initialize tot_steps_, depending on total number of steps and border_weight
  tot_steps_=steps_beta_*steps_pres_+(border_weight_-1)*2*(steps_beta_+steps_pres_);
  tot_steps_*=tot_umbrellas_;
  tot_steps_+=P0_contribution_; //FIXME is this the best way to add it?

//print some info
  log.printf("  Total temp steps (in beta) = %d\n",steps_beta_);
  for(unsigned i=0; i<steps_beta_; i++)
    log.printf("   %2d. beta=% .10f  temp=% g\n",i+1,beta_p_[i],1./(plumed.getAtoms().getKBoltzmann()*beta_p_[i]));
  if(NumParallel_>steps_beta_)
    log.printf(" +++ WARNING +++ number of parallel threads is greater than beta steps. Using SERIAL might be faster\n");
  log.printf("  Total pres steps = %d\n",steps_pres_);
  for(unsigned j=0; j<steps_pres_; j++)
    log.printf("   %2d. pres=% .10f\n",j+1,pres_p_[j]);
  //this are actually set from the beginning
  log.printf("  Total number of umbrellas = %u\n",tot_umbrellas_);
  log.printf("    with SIGMA = %g\n",sigma_);
  log.printf("    in CV range [%g,%g]\n",center_[0],center_[tot_umbrellas_-1]);
  log.printf("  Total deltaFs = %u\n",steps_beta_*steps_pres_*tot_umbrellas_);
  if(border_weight_!=1)
    log.printf(" --- using a weight different from 1 for the temperature-pressure border, BORDER_WEIGHT = %g ---\n",border_weight_);
}

void OPESmultiPhase::init_from_obs()
{
//in case of multiple walkers gather all the statistics
  if(NumWalkers_>1)
  {
    obs_steps_*=NumWalkers_;
    std::vector<double> all_obs_ene(obs_steps_);
    std::vector<double> all_obs_vol(obs_steps_);
    std::vector<double> all_obs_cv(obs_steps_);
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Allgather(obs_ene_,all_obs_ene);
      multi_sim_comm.Allgather(obs_vol_,all_obs_vol);
      multi_sim_comm.Allgather(obs_cv_,all_obs_cv);
    }
    comm.Bcast(all_obs_ene,0);
    comm.Bcast(all_obs_vol,0);
    comm.Bcast(all_obs_cv,0);
    obs_ene_=all_obs_ene;
    obs_vol_=all_obs_vol;
    obs_cv_=all_obs_cv;
  }

//estimate steps from Neff and initialize integration grid
  if(steps_beta_==0) //beta_p_[1]=max_beta->min_temp and beta_p_[0]=min_beta->max_temp
    steps_beta_=estimate_steps(beta_p_[1]-beta_,beta_p_[0]-beta_,obs_ene_,"TEMP");
  if(steps_pres_==0) //pres_p_[0]=min_pres, pres_p_[1]=max_pres
    steps_pres_=estimate_steps(beta_*(pres_p_[0]-pres_),beta_*(pres_p_[1]-pres_),obs_vol_,"PRESSURE");
  init_integration_grid();

//initialize deltaF_
  deltaF_.resize(steps_beta_,std::vector< std::vector<double> >(steps_pres_,std::vector<double>(tot_umbrellas_)));
  std::vector<double> first_umbrellas(tot_umbrellas_);
  for(unsigned k=0; k<tot_umbrellas_; k++)
  {
    first_umbrellas[k]=0.5*std::pow((obs_cv_[0]-center_[k])/sigma_,2)/beta_;
    if(first_umbrellas[k]>barrier_)
      first_umbrellas[k]=barrier_;
  }
  for(unsigned i=0; i<steps_beta_; i++)
    for(unsigned j=0; j<steps_pres_; j++)
      for(unsigned k=0; k<tot_umbrellas_; k++)
        deltaF_[i][j][k]=(beta_p_[i]/beta_-1.)*obs_ene_[0]+(beta_p_[i]/beta_*pres_p_[j]-pres_)*obs_vol_[0]+first_umbrellas[k];
  for(unsigned i=0; i<steps_beta_; i++)
    for(unsigned j=0; j<steps_pres_; j++)
      for(unsigned k=0; k<tot_umbrellas_; k++)
        for(unsigned t=1; t<obs_steps_; t++) //starts from t=1
          deltaF_[i][j][k]+=-1./beta_*std::log1p(get_weight(beta_p_[i],obs_ene_[t],pres_p_[j],obs_vol_[t],center_[k],obs_cv_[t],deltaF_[i][j][k])/t)-1./beta_*std::log1p(-1./(1.+t));
  if(calc_work_)
    old_deltaF_=deltaF_;
  obs_ene_.clear();
  obs_vol_.clear();
  obs_cv_.clear();

//print initialization to file
  deltaFsOfile_.printField("time",getTime());
  deltaFsOfile_.printField("rct",rct_);
  for(unsigned i=0; i<steps_beta_; i++)
    for(unsigned j=0; j<steps_pres_; j++)
      for(unsigned k=0; k<tot_umbrellas_; k++)
        deltaFsOfile_.printField("deltaF_"+std::to_string(i+1)+"_"+std::to_string(j+1)+"_"+std::to_string(k+1),deltaF_[i][j][k]);
  deltaFsOfile_.printField();
}

unsigned OPESmultiPhase::estimate_steps(const double left_side,const double right_side,const std::vector<double>& obs,const std::string msg) const
{ //uses Neff to estimate the grid spacing
  if(left_side==0 && right_side==0)
  {
    log.printf(" +++ WARNING +++ MIN_%s and MAX_%s are equal to %s, using single step\n",msg.c_str(),msg.c_str(),msg.c_str());
    return 1;
  }
  auto get_neff_HWHM=[](const double side,const std::vector<double>& obs,const double av_obs) //HWHM = half width at half maximum. neff is in general not symmetric
  {
    //func: Neff/N-0.5 is a function between -0.5 and 0.5
    auto func=[](const long double delta,const std::vector<double> obs, const double av_obs)
    {
      long double sum_w=0;
      long double sum_w2=0;
      for(unsigned t=0; t<obs.size(); t++)
      {
        const long double w=std::exp(-delta*(obs[t]-av_obs));
        sum_w+=w;
        sum_w2+=w*w;
      }
      return sum_w*sum_w/sum_w2/obs.size()-0.5;
    };
    //here we find the root of func using the regula falsi (false position) method
    //but any method would be OK, not much precision is needed. src/tools/RootFindingBase.h looked complicated
    const double tolerance=1e-4; //seems to be a good default
    double a=0; //default is right side case
    double func_a=0.5;
    double b=side;
    double func_b=func(side,obs,av_obs);
    if(func_b>=0)
      return 0.0; //no zero is present!
    if(b<0)//left side case
    {
      std::swap(a,b);
      std::swap(func_a,func_b);
    }
    double c=a;
    double func_c=func_a;
    while(std::abs(func_c)>tolerance)
    {
      if(func_a*func_c>0)
      {
        a=c;
        func_a=func_c;
      }
      else
      {
        b=c;
        func_b=func_c;
      }
      c=(a*func_b-b*func_a)/(func_b-func_a);
      func_c=func(c,obs,av_obs);//func is evaluated only here, it might be expensive
    }
    return std::abs(c);
  };

//set average to zero, for numerical stability
  double av_obs=0;
  for(unsigned t=0; t<obs.size(); t++)
    av_obs+=obs[t];
  av_obs/=obs.size();

//estimation
  double left_HWHM=0;
  if(left_side!=0)
    left_HWHM=get_neff_HWHM(left_side,obs,av_obs);
  double right_HWHM=0;
  if(right_side!=0)
    right_HWHM=get_neff_HWHM(right_side,obs,av_obs);
  if(left_HWHM==0)
  {
    right_HWHM*=2;
    if(left_side==0)
      log.printf(" --- MIN_%s is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ MIN_%s is very close to %s\n",msg.c_str(),msg.c_str());
  }
  if(right_HWHM==0)
  {
    left_HWHM*=2;
    if(right_side==0)
      log.printf(" --- MAX_%s is equal to %s\n",msg.c_str(),msg.c_str());
    else
      log.printf(" +++ WARNING +++ MAX_%s is very close to %s\n",msg.c_str(),msg.c_str());
  }
  if(left_HWHM==0 && right_HWHM==0)
  {
    log.printf(" +++ WARNING +++ %s range is very narrow, using MIN_%s and MAX_%s as only steps\n",msg.c_str(),msg.c_str(),msg.c_str());
    return 2;
  }
  const double grid_spacing=left_HWHM+right_HWHM;
  log.printf("  Estimated %s spacing (with beta) = %g\n",msg.c_str(),grid_spacing);
  unsigned steps=std::ceil(std::abs(right_side-left_side)/grid_spacing);
  plumed_massert(steps>1,"something went wrong and estimated grid spacing for "+msg+" gives a step="+std::to_string(steps));
  return steps;
}

}
}
