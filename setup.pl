#!/usr/bin/perl

# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2014-2016

# This is a pre-processor used to generate the fortran code for your chosen shell model. It generates shellmodels.f90 with your model setup using the template shellmodels_plain.f90.

use warnings;

#-----------------------------------
# MODEL TYPES
#-----------------------------------
# Model 1-4    Helical Sabra models    (unpublished master thesis: N. M. Rathmann (2014), De Pietro et al., 2015) 
# Model 11-14  New helical models      (Rathmann and Ditlevsen, 2016), referred to as the RD model below.
# Model 15     New mixed helical model (Rathmann and Ditlevsen, 2016), referred to as the RD model below.
#-----------------------------------

my $MODEL        = $ARGV[0];
my $LAMBDA_m     = $ARGV[1]; # Lambda used to determine epsilon and xi coeficients.
my $LAMBDA_g     = $ARGV[2]; # Lambda used to determine helical interaction weight.
my $LAMBDA_G     = $ARGV[3]; # Lambda used to determine shape weight.
my $DEMARC_MODEL = $ARGV[4]; # Helical sub-interaction to demarcate.
my $DEMARC_VALUE = $ARGV[5]; # Helical sub-interaction bias towards model $DEMARC_MODEL: 0 produces bias, 1 turns of all sub-interactions other than $DEMARC_MODEL.
my $FIXED_P      = $ARGV[6]; 
my $FIXED_Q      = $ARGV[7];
my $NSH          = $ARGV[8]; 
my $F_IN         = $ARGV[9];
my $F_OUT        = $ARGV[10];

print "DEBUG: lambda_m = $LAMBDA_m, lambda_g = $LAMBDA_g, lambda_G = $LAMBDA_G\n";

#--------------------------------------------------------------------------

# Model 1-4 Sp and Sq values (s^prime and s^double-prime)
my @Sp = (undef, -1, -1, +1, +1); # Set first zero'th entry to undef so that indexes match the submodel numbering.
my @Sq = (undef, +1, -1, -1, +1);

# Model triple-correlator parity
my $SIGN_P = '+';
my $SIGN_N = '-';
if ($MODEL < 10) { $SIGN_N = '+'; }   
my $IS_SABRA = ($SIGN_N eq '+');

# Submodels simulated
my @MODELLIST;
if    ($MODEL == 15) {@MODELLIST = 1..4;}        # New model, all submodels coupled
elsif ($MODEL  > 10) {@MODELLIST = ($MODEL-10);} # New model, single submodels
elsif ($MODEL  >  0) {@MODELLIST = ($MODEL);}    # Helical Sabra submodels, uncoupled
print sprintf("Number of helical interaction types (submodels) simulated: %i\n\n", scalar(@MODELLIST));

# Some flags
my $IS_MULTI_SUBMODEL = (scalar(@MODELLIST) > 1);
my $IS_PQ_FIXED = !(($FIXED_P == 0) && ($FIXED_Q == 0));

# Weights
sub G { my ($L,$p,$q) = @_;         return (1/8)*( 2*$L**(-2*$q) + 2*$L**(-2*$p) + 2 - $L**(-2*($p+$q)) - $L**(2*($p-$q)) - $L**(2*($q-$p)) )**0.5; }
sub g { my ($L,$p,$q,$Sp,$Sq) = @_; return -1*$Sp*$Sq*(1 + $Sp*$L**$p - $Sq*$L**$q)*($Sp*$L**$p - $Sq*$L**$q); }
#sub g { my ($L,$p,$q,$Sp,$Sq) = @_; return abs(-1*$Sp*$Sq*(1 + $Sp*$L**$p - $Sq*$L**$q)*($Sp*$L**$p - $Sq*$L**$q)); }

sub eps_RD    { my ($L,$p,$q,$Sp,$Sq) = @_; return ($Sp - $Sp*$Sq*$L**($q))/($Sp*$L**($p) - $Sq*$L**($q)); }
sub xi_RD     { my ($L,$p,$q,$Sp,$Sq) = @_; return $Sp*$Sq*eps_RD($L,$p,$q,$Sp,$Sq) - $Sq; }
sub eps_Sabra { my ($L,$p,$q,$Sp,$Sq) = @_; return $Sp*eps_RD($L,$p,$q,$Sp,$Sq); }
sub xi_Sabra  { my ($L,$p,$q,$Sp,$Sq) = @_; return 1 - eps_Sabra($L,$p,$q,$Sp,$Sq); }
my $intactcoefs = { eps => !$IS_SABRA ? \&eps_RD : \&eps_Sabra , xi => !$IS_SABRA ? \&xi_RD : \&xi_Sabra };

#--------------------------------------------------------------------------
# DETERMINE INTERACTION SCOPE (maximum possible q-value for valid geometry)
#--------------------------------------------------------------------------

sub test_triineq    { my ($L,$p,$q) = @_; return (1 + $L**($p) >= $L**($q)); }
#sub test_myrestrict { my ($p,$q,$FIXED_P,$FIXED_Q) = @_; return (!$IS_PQ_FIXED) || (($p == $FIXED_P) && ($q == $FIXED_Q));} # Pick out only specific interactions.
sub test_myrestrict { 
        my ($p,$q,$FIXED_P,$FIXED_Q) = @_; 
        my $FOUND_SPECIFIED_GEOM = 0;
        my @Plist = split(',',$FIXED_P);
        my @Qlist = split(',',$FIXED_Q);
        foreach my $ii (0..$#Plist) {
                if ( ($Plist[$ii] eq $p) && ($Qlist[$ii] eq $q) ) {$FOUND_SPECIFIED_GEOM = 1;} # Pick out only specified interactions.
        }
        return (!$IS_PQ_FIXED) || $FOUND_SPECIFIED_GEOM; 
} 

my ($q_max, $num_of_triads) = (0,0);
my ($G_sum, $g_sum, $Gg_sum) = (0,0,0); # For normalisation
my @plist = (); my @qlist = ();

print "Allowed geometries:\n";
for (my $q=2; $q <= 1000; $q++) { # Use some random large upper bound
for (my $p=1; $p <= $q-1; $p++) {

        if ( test_triineq($LAMBDA_G,$p,$q) && test_myrestrict($p,$q,$FIXED_P,$FIXED_Q) ) {
        
                my $G = G($LAMBDA_G,$p,$q);
                $G_sum += $G;
               
                print "    p,q=$p,$q   G=".sprintf('%.2e',$G)."\n";
               
                foreach my $submodel (@MODELLIST) { 
                        my $g = g($LAMBDA_g,$p,$q,$Sp[$submodel],$Sq[$submodel]);
                        $g_sum  += $g; 
                        $Gg_sum += $G*$g; 
                }
                
                $q_max = $q;
                $num_of_triads++;
                push(@plist,$p); push(@qlist,$q);
        }
}
}

my $IS_MULTI_GEOM = ($q_max > 2) && !$IS_PQ_FIXED;
if (($num_of_triads>1) && !$IS_MULTI_GEOM) {print "\nERROR: More triads found that expected => Aborting.\n"; die();}
print sprintf('Weight sums: sum(G) = %.3f, sum(g) = %.3f, sum(G*g) = %.3f',$G_sum,$g_sum,$Gg_sum)."\n\n";

#--------------------------------------------------------------------------
# DETERMINE REPLACEMENTS THIS PRE-PROCESSOR SHOULD MAKE
#--------------------------------------------------------------------------

# These are the array (parameter) constructors.
my (@G, @g, @eps, @xi, @d1, @d2, @d3);
my ($INTACT_P, $INTACT_N) = ('','');
my @CORR_P = ('','','','');my @CORR_N = ('','','','');my @ECORR  = ('','','','');
my ($test_G_sum, $test_g_sum, $test_Gg_sum) = (0,0,0);

print sprintf("Using triads (normalising G -> G/%.2e):\n",$Gg_sum);
foreach my $triad (1..$num_of_triads) {

        my $p = $plist[$triad-1];
        my $q = $qlist[$triad-1];
        if ($p>1 && $q>1) { print "\n"; } # For nicer visual output
        
        #--------------------------------------------------------------------------
        # d(u_n)/dt INTERACTIONS
        #--------------------------------------------------------------------------
        $INTACT_P .= " $SIGN_P".(($IS_SABRA) ? 'cmplx(0,1)*' : '').'('.  join('+', map { "Gg($triad,$_)*(DUDT_P_$_(STEP,N,$triad,$p,$q))" } @MODELLIST ) .')';
        $INTACT_N .= " $SIGN_N".(($IS_SABRA) ? 'cmplx(0,1)*' : '').'('.  join('+', map { "Gg($triad,$_)*(DUDT_N_$_(STEP,N,$triad,$p,$q))" } @MODELLIST ) .')';
        
        #--------------------------------------------------------------------------
        # WEIGHTS AND COEFICIENT MATRICES
        #--------------------------------------------------------------------------
        $G[$triad] = G($LAMBDA_G,$p,$q)/$Gg_sum; # Normalise by $Gg_sum
        $test_G_sum += $G[-1];

        foreach my $submodel (1..4) {

                my $IS_SUBMODEL_ENABLED = scalar( grep { $submodel eq $_ } @MODELLIST );
        
                $eps[$triad][$submodel] = $intactcoefs->{'eps'}($LAMBDA_m,$p,$q,$Sp[$submodel],$Sq[$submodel]);
                $xi[$triad][$submodel]  = $intactcoefs->{'xi'}( $LAMBDA_m,$p,$q,$Sp[$submodel],$Sq[$submodel]);
                my $_g   = '.0';
                
                if ($IS_SUBMODEL_ENABLED) {
                
                        # Regular demarcation - turns OFF the three interaction types which are NOT $DEMARC_MODEL
                        if ($DEMARC_VALUE >= 0) { 
                        
                                $_g = g($LAMBDA_g,$p,$q,$Sp[$submodel],$Sq[$submodel]);
                                if ($submodel == $DEMARC_MODEL) {
                                        foreach my $j (1..4) { if ($j != $submodel) { $_g += $DEMARC_VALUE*g($LAMBDA_g,$p,$q,$Sp[$j],$Sq[$j]); } }
                                } 
                                else { $_g *= (1-$DEMARC_VALUE); }
                        }
                        
                        # Reverse demarcation - turns OFF the specific interaction type
                        elsif ($DEMARC_VALUE < 0) { 
                        
                                $_g = g($LAMBDA_g,$p,$q,$Sp[$submodel],$Sq[$submodel]);
                                if ($submodel == $DEMARC_MODEL) { $_g *= (1-abs($DEMARC_VALUE)); } 
                                else {                     $_g += abs($DEMARC_VALUE/3)*g($LAMBDA_g,$p,$q,$Sp[$DEMARC_MODEL],$Sq[$DEMARC_MODEL]); }
                        }
                }
                
                $g[$triad][$submodel] = $_g;
                $test_g_sum  += $_g;
                $test_Gg_sum += $G[-1]*$_g;

                print "    p,q=$p,$q (".(($Sp[$submodel]>0) ? '+' : '-').",".(($Sq[$submodel]>0) ? '+' : '-')."):" . "\tg=".(($_g != 0) ? sprintf('%2.2e',$_g) : 0). ", G=".(($G[-1] != 0 ) ? sprintf('%.2e',$G[-1]) : 0). "\tG*g=". sprintf('%2.2e',$G[-1]*$_g). "\n";
                
                my $zeros_qmax = "(.0,ii=1,$q_max)"; # Set interaction coefs of zero-velocity boundary shells to zero (to be safe).
                my $zeros_full = "(.0,ii=1,".($NSH+2*$q_max).")"; # Set interaction coefs of zero-velocity boundary shells to zero (to be safe).
                my $n0   = $q_max+1;
                my $nend = $q_max+$NSH;
                if ($IS_SUBMODEL_ENABLED) {
                        $d1[$triad][$submodel] = "$zeros_qmax,(k(ii),ii=$n0,$nend),$zeros_qmax";
                        $d2[$triad][$submodel] = "$zeros_qmax,(-k(ii)*".$eps[$triad][$submodel]."/($LAMBDA_m**$p),ii=$n0,$nend),$zeros_qmax";
                        $d3[$triad][$submodel] = "$zeros_qmax,(k(ii)*".$xi[$triad][$submodel]."/($LAMBDA_m**$q),ii=$n0,$nend),$zeros_qmax";
                } else {
                        $d1[$triad][$submodel] = $zeros_full;
                        $d2[$triad][$submodel] = $zeros_full;
                        $d3[$triad][$submodel] = $zeros_full;
                }
        }
        
        #--------------------------------------------------------------------------
        # CORRELATOR AGGREGATIONS
        #--------------------------------------------------------------------------
        foreach my $submodel (@MODELLIST) {
                $CORR_P[$submodel-1] .= "corr_p(nn,$triad,$submodel)=corr_p(nn,$triad,$submodel)+CORR_P_$submodel(nn,$p,$q); ";
                $CORR_N[$submodel-1] .= "corr_n(nn,$triad,$submodel)=corr_n(nn,$triad,$submodel)+CORR_N_$submodel(nn,$p,$q); ";
                $ECORR[$submodel-1]  .= "Ecorr(nn,$triad,$submodel)=CORR_P_$submodel(nn,$p,$q)+MODEL_PARITY*CORR_N_$submodel(nn,$p,$q); ";
        }
}

print sprintf("Weight sums sanity check: sum(G) = %.3f, sum(g) = %.3f, sum(G*g) = %.3f",$test_G_sum,$test_g_sum,$test_Gg_sum);
my $precision = 0.000001;
if ($test_Gg_sum < 1+$precision && $test_Gg_sum > 1-$precision ) {print " ... OK (sum(G*g) = 1)\n"}
else {die('*** NORMALISATION FAILED')}

#--------------------------------------------------------------------------

my @G_joinme = (); my @g_joinme = (); my @eps_joinme = (); my @xi_joinme = (); my @d1_joinme = (); my @d2_joinme = (); my @d3_joinme = ();
foreach my $triad (1..$num_of_triads) {
        push(@G_joinme,   $G[$triad]);
}
foreach my $submodel (1..4) {
foreach my $triad (1..$num_of_triads) {
        push(@g_joinme,   sprintf('%.8f',$g[$triad][$submodel]));
        push(@eps_joinme, $eps[$triad][$submodel]);
        push(@xi_joinme,  $xi[$triad][$submodel]);
        push(@d1_joinme,  $d1[$triad][$submodel]);
        push(@d2_joinme,  $d2[$triad][$submodel]);
        push(@d3_joinme,  $d3[$triad][$submodel]);
}}

my %REPLACE;
$REPLACE{'REPLACE__MODEL'}  = $MODEL;
$REPLACE{'REPLACE__NSH'}    = $NSH;
$REPLACE{'REPLACE__LAMBDA'} = sprintf('%.8f',$LAMBDA_m); # Must not be integer, convert to float
$REPLACE{'REPLACE__q_MAX'}  = $q_max;
$REPLACE{'REPLACE__NUM_TRIAD_GEOMS'} = $num_of_triads;
$REPLACE{'REPLACE__SUBMODEL_CLASS'} = ($IS_MULTI_SUBMODEL) ? 5 : $MODELLIST[0];
$REPLACE{'REPLACE__DEMARCATION_MODEL'} = $DEMARC_MODEL;
$REPLACE{'REPLACE__DEMARCATION_VALUE'} = $DEMARC_VALUE;
$REPLACE{'REPLACE__DUDT_P'} = $INTACT_P;
$REPLACE{'REPLACE__DUDT_N'} = $INTACT_N;
$REPLACE{'REPLACE__p_LIST'} = '['.join(',',@plist).']';
$REPLACE{'REPLACE__q_LIST'} = '['.join(',',@qlist).']';
$REPLACE{'REPLACE__G_CONSTRUCTOR'}   = '['.join(',',@G_joinme)."]"; 
$REPLACE{"REPLACE__g_CONSTRUCTOR"}   = '['.join(',',@g_joinme).']';
$REPLACE{"REPLACE__eps_CONSTRUCTOR"} = '['.join(',',@eps_joinme).']';
$REPLACE{"REPLACE__xi_CONSTRUCTOR"}  = '['.join(',',@xi_joinme).']';
$REPLACE{"REPLACE__d1_CONSTRUCTOR"}  = '['.join(',',@d1_joinme).']';
$REPLACE{"REPLACE__d2_CONSTRUCTOR"}  = '['.join(',',@d2_joinme).']';
$REPLACE{"REPLACE__d3_CONSTRUCTOR"}  = '['.join(',',@d3_joinme).']';

foreach my $submodel (1..4) {
        $REPLACE{"REPLACE__SUM_CORR_P_$submodel"} = $CORR_P[$submodel-1]."\n";
        $REPLACE{"REPLACE__SUM_CORR_N_$submodel"} = $CORR_N[$submodel-1]."\n";
        $REPLACE{"REPLACE__ECORR_$submodel"}      = $ECORR[$submodel-1]."\n";
}

#--------------------------------------------------------------------------
# SEARCH AND REPLACE
#--------------------------------------------------------------------------

open(my $H_IN,  '<', $F_IN)  or die "Could not open file '$F_IN' $!";
open(my $H_OUT, '>', $F_OUT) or die "Could not open file '$F_OUT' $!";
while (my $row = <$H_IN>) {
        foreach my $key (keys %REPLACE) {$row =~ s/$key/$REPLACE{$key}/;}
        print $H_OUT $row;
}
print "\n";
