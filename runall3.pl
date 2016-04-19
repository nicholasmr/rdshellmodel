#!/usr/bin/perl
# Nicholas Mossor Rathmann, 2014-2015

# SIMULATION TIME FOR STATS CONVERGENCE IN ALL MODELS
#       our $NT     = 1000*2000;
#       our $NTI    = 1000*1000; (with dt=1e-8) # Inner loop should be 1000*1000 for dt = 1e-8, or 1000*100 for dt = 1e-7 (gives accurate sol. also)

use warnings;
use strict;

#------------ DEBUG ----------------

our $DEBUG__NO_MAKE   = 0;
our $DEBUG__NO_SUBMIT = 0;

#------------ GENERAL SETUP --------

our $MULTITHREADED = 0;
our $CORES         = 6;

#our $DT  = '1.0D-7';
#our $NTI = 1000*100; # use ~ 1000*1000 with $DT  = '1.0D-8'

our $DT  = '7.0D-8';
our $NTI = 1000*1000/7; # use ~ 1000*1000 with $DT  = '1.0D-8'

our $FMAGNITUDE = 1;

our @VISC_LARGESCALE;
our @VISC_SMALLSCALE;

$VISC_LARGESCALE[12] = "1.0D3"; # for k^-4
#$VISC_LARGESCALE[12] = "1.0D2"; # for k^-2
$VISC_LARGESCALE[11] = $VISC_LARGESCALE[12]; # for k^-2
$VISC_LARGESCALE[13] = $VISC_LARGESCALE[12]; # for k^-2
$VISC_LARGESCALE[14] = $VISC_LARGESCALE[12]; # for k^-2

$VISC_SMALLSCALE[12] = '1.0D-11';
$VISC_SMALLSCALE[11] = $VISC_SMALLSCALE[12];
$VISC_SMALLSCALE[13] = $VISC_SMALLSCALE[12];
$VISC_SMALLSCALE[14] = $VISC_SMALLSCALE[12];

#----------------------------
our ($LAMBDA_G, $LAMBDA_m, $NSH, $DTYPE, $FTYPE, $FSH, $DFSH, $FSH_MID, $FSH_TOP, $NT, $DEMARC_MODEL, $DEMARC_VALUE, $FIXED_P, $FIXED_Q, $DO_AGGREGATIONS, $RESUME_AGGREGATIONS, $RESUME, $CPRESUME, $CPINIT, $LBL, $BASEDIR);

my (@plist,@qlist);

$LAMBDA_G = 1.10;
$LAMBDA_m = '1.10D0';
#@plist = (1,6,11,15,21);
#@qlist = (2,7,12,16,22);

@plist = (1,11,21);
@qlist = (2,12,22);

$FTYPE    = +2;
#$FTYPE    = -2;
$DTYPE    = 3;

my $NSH_ref = 143 + 2*40;
my $FSH_ref = 48  + 1*30;

$CPINIT = 0;
$DEMARC_VALUE = 0; $DEMARC_MODEL = 0;
$FIXED_P = 0; $FIXED_Q = 0;
$LBL = 'paper3'; my $log_a = sprintf("NEWEXP");


#-----------------------------------

$RESUME = 0;

if ($RESUME == 1) {$NT = 1000*11;}
else              {$NT = 1000*10;}
$DO_AGGREGATIONS = $RESUME;
$CPRESUME = $RESUME;
$RESUME_AGGREGATIONS = 0;

foreach my $i (0..$#qlist) {


        $FIXED_P = $plist[$i];
        $FIXED_Q = $qlist[$i];
        $DFSH = $FIXED_P; 

        $BASEDIR  = sprintf("../../paper3/pq%i%i",$FIXED_P,$FIXED_Q);

        # Very high Re
        $NSH = $NSH_ref; $FSH = $FSH_ref; 
        launch(11,  $FTYPE, $log_a);                       
        launch(12,  $FTYPE, $log_a);
        launch(13,  $FTYPE, $log_a); 
        launch(14,  $FTYPE, $log_a); 
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

sub launch {

        my ($m, $ftype, $label) = @_;
        
        my $lambda = $LAMBDA_m; $lambda =~ s/(\.0|D0)//g; eval('$lambda = '.$lambda); 
	my $FEXT   = sprintf('m%i_nsh%i_lam%1.1f_dtype%i_ftype%i_fsh%i',$m,$NSH,$lambda, $DTYPE,$ftype,$FSH);
        my $OUTDIR = "$BASEDIR/$FEXT";

        ################
        # Makefile stuff
        ################
        my $SETUP_CONFIG = "MODEL=\"$m\" NSH=\"$NSH\" LAMBDA_m=\"$lambda\" LAMBDA_g=\"$lambda\" LAMBDA_G=\"$LAMBDA_G\" DEMARC_MODEL=\"$DEMARC_MODEL\" DEMARC_VALUE=\"$DEMARC_VALUE\" FIXED_P=\"$FIXED_P\" FIXED_Q=\"$FIXED_Q\"";
        my $OTHER = "OUTDIR=\"$OUTDIR\" MULTITHREADED=\"$MULTITHREADED\" RESUME=\"$RESUME\" DO_AGGREGATIONS=\"$DO_AGGREGATIONS\" RESUME_AGGREGATIONS=\"$RESUME_AGGREGATIONS\"";
        #----------------
        my $MACROS__MAIN__FORCE    = "-DFTYPE=$ftype -DFSH=$FSH -DDFSH=$DFSH -DFMAGNITUDE=$FMAGNITUDE";
        my $MACROS__MAIN__DISSP    = "-DDTYPE=$DTYPE -DVISC_SMALLSCALE=$VISC_SMALLSCALE[$m] -DVISC_LARGESCALE=$VISC_LARGESCALE[$m]";
        my $MACROS__MAIN__NUMERICS = "-DNT=$NT -DNTI=$NTI -DDT=$DT";
        #----------------
        my $MAKEFILE_ARGS = "$OTHER $SETUP_CONFIG MACROS__MAIN__FORCE=\"$MACROS__MAIN__FORCE\" MACROS__MAIN__DISSP=\"$MACROS__MAIN__DISSP\"  MACROS__MAIN__NUMERICS=\"$MACROS__MAIN__NUMERICS\" ";
	my $CMD = "make shellmodel $MAKEFILE_ARGS EXT=\"__$FEXT\" ";
	if (!$DEBUG__NO_MAKE) { print $CMD; system($CMD); }
        #----------------
        
        ################
        # SLURM stuff
        ################
        
        if ($CPRESUME) { system("mv $OUTDIR/dump__$FEXT.nc   $OUTDIR/resume__$FEXT.nc") }
        if ($CPINIT)   { system("mv $OUTDIR/resume__$FEXT.nc $OUTDIR/init__$FEXT.nc")   }

	my $THENODE = "-w node04"; $THENODE = "";
        
	my $SBATCHFILE = "${OUTDIR}/$FEXT.sh";
	open(JOB, ">$SBATCHFILE");
	print JOB "#!/bin/bash			\n";
	print JOB "#SBATCH -p flex	        \n";

	if ($MULTITHREADED) {
	        # this option is needed for multithreaded (e.g. OpenMP) jobs, it tells SLURM to allocate N cores per task allocated; typically N should be equal to the number of threads you program spawns, e.g. it should be set to the same number as OMP_NUM_THREADS
#                print JOB "#SBATCH --ntasks-per-core=1  \n";
                print JOB "#SBATCH --cpus-per-task=$CORES \n";	 
	        print JOB "export OMP_NUM_THREADS=$CORES\n";
	} else {
        	print JOB "#SBATCH -s \n";
	}
	print JOB "cd ${OUTDIR}\n";
	print JOB "srun -s -u -o $FEXT.out shellmodel__$FEXT\n";
	close(JOB);

        my $SBATCH = "sbatch $THENODE $SBATCHFILE"; print "$SBATCH\n";
	if (!$DEBUG__NO_SUBMIT) { system($SBATCH); print "\n--------------------------\n"; }
	else {print "*******************\n NOT SUBMITTED TO SLURM \n*******************\n";}
	
        open(STATUSLOG,  '>>', "STATUS_${label}.sh");
        print STATUSLOG "perl printstatus.pl ".sprintf(' %s',$OUTDIR).";\n";
        close(STATUSLOG);
}

sub clearlog {
        my ($label) = @_;
        open(STATUSLOG, '>', "STATUS_${label}.sh"); close(STATUSLOG);
}
