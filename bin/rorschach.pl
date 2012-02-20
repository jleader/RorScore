#!perl

use strict;
use warnings;

# Example input (tab-separated, lf or cr/lf delimited lines?)
# Card \t Resp # \t Location and DQ \t Loc # \t Determinant(s) and FQ \t
#     (2) \t Content(s) \t Pop \t Z-score \t Spec. Scores
# I \t 1 \t W+ \t  \t F- \t 2 \t Fd \t  \t 4 \t "DV, MOR"
# \t 2 \t Wo \t  \t Fo \t  \t Ad \t  \t 1 \t MOR
# (spaces added around the tabs for easier readability)

print main(shift) unless caller;

sub main {
    my ($infname) = @_;
    my @romans = qw(0 I II III IV V VI VII VIII IX X);
    my %romans;
    foreach my $i (0..$#romans) {
        $romans{$romans[$i]} = $i;
    }

    my @determinants = (qw(M FM m FC CF C Cn), "FC'", "C'F", "C'",
        qw(FT TF T FV VF V FY YF Y Fr rF FD F));
    my %isdet = map { $_ => 1 } @determinants;

    my @contents = ('H', '(H)', 'Hd', '(Hd)', 'Hx', 'A', '(A)', 'Ad', '(Ad)',
        qw(An Art Ay Bl Bt Cg Cl Ex Fd Fi Ge Hh Ls Na Sc Sx Xy Id));
    my %contents = map { $_ => 1 } @contents;

    my %special1 = (DV=>1, INC=>2, DR=>3, FAB=>4, ALOG=>5, CON=>7);
    my @special1 = keys %special1;
    my %special2 = (DV=>2, INC=>4, DR=>6, FAB=>7);
    my @special0 = qw(AB AG COP CP GHR PHR MOR PER PSV);
    my %special = (%special1,
                   (map { $_ . '2' => $special2{$_} } keys %special2),
                   (map { $_ => 0 } @special0));
    my %specialsynonyms = (
        CONTAM  => 'CON',
        INCOM1  => 'INC',
        INCOM2  => 'INC2',
        FABCOM1 => 'FAB',
        FABCOM2 => 'FAB2',
        map { $_ . '1' => $_ } @special1,
    );

    my ($hasresp, $cardnum, $respnum);
    my @resps = (0) x 11;
    my ($Zf, $ZSum, $W, $D, $Dd, $S, %DQ, %FQ, %MQ, %WDQ,
        %ap, %Map, @blends, $twos,
        @approaches,
        $raw6, $weight6, $popular);
    my %ccounts = map { $_ => 0 } @contents;
    my %sc = map { $_ => 0 }
             (@special1, @special0, map { $_ . '2' } keys %special2);
    my $sminus = 0;
    my %singles = map { $_ => 0 } @determinants;
    my %mults = map { $_ => 0 } @determinants;
    open my $inf, '<', $infname or die "Unable to read from $infname: $!";
    while (<$inf>) {
        my @fields = split /\t/;
        #print "DEBUG: fields=<@fields>\n";
        if ($. == 1 && $fields[0] =~ /^ \s* card /ix) {
            $hasresp = !!($fields[1] =~ /^ \s* resp /ix);
            #print "DEBUG hasresp=$hasresp\n";
            next;
        }
        else {
            ++$respnum;
            if ($fields[0] !~ /^ \s* $/x) {
                # new card
                #print "DEBUG: field0=$fields[0]\n";
                my $newcardnum = $romans{$fields[0]} || $fields[0];
                ++$cardnum;
                die "Mis-orderd cards: $fields[0], expected $romans[$cardnum]"
                    unless $cardnum == $newcardnum;
                $cardnum = $newcardnum;
                #print "DEBUG: new card $cardnum\n";
            }
            ++$resps[$cardnum];
            die 'Missing first card number' unless $cardnum;
            die "Bad response number $fields[1], expected $respnum"
                if $hasresp && $respnum != $fields[1];
            splice @fields, 0, ($hasresp ? 2 : 1);

            # strip whitespace and enclosing double-quotes
            @fields = map { s/\s+//g; s/^"(.*)"$/$1/; $_ } @fields;
            #print "DEBUG: clean fields=<@fields>\n";

            my ($location, $locno, $detsfq, $two, $content, $pop,
                $zscore, $specscores, @extras) = @fields;
            warn scalar(@extras) . ' extra fields' if @extras;

            $location =~ m{^ (W|D|Dd) (S?) (\+|o|v/\+|v) $}x
                or die 'Bad location/DQ column';
            my ($loc, $loc_s, $dq) = ($1, $2, $3);

            my $hasW  = !!($loc eq 'W');
            my $hasD  = !!($loc eq 'D');
            my $hasDd = !!($loc eq 'Dd');
            my $hasS  = !!($loc_s eq 'S');
            $W  += $hasW;
            $D  += $hasD;
            $Dd += $hasDd;
            $S  += $hasS;

            $DQ{p}  += !!($dq eq '+');
            $DQ{o}  += !!($dq eq 'o');
            $DQ{vp} += !!($dq eq 'v/+');
            $DQ{v}  += !!($dq eq 'v');

            push @{$approaches[$cardnum]}, $loc . $loc_s;

            $detsfq =~ /^ ([aCDFmMn'prTVY.]*) ([+ou-]?) $/x
                or die 'Bad determinants/FQ column';
            my ($dets, $fq) = ($1, $2);
            my @dets = split /\./, $dets;
            foreach (@dets) {
                # strip out and count superscript a or p
                if (s/^ (m|M|FM) ([ap]) $/$1/x) {
                    ++$ap{$2};
                    ++$Map{$2} if $1 eq 'M';
                }
            }
            my @baddets = grep { !$isdet{$_} } @dets;
            die 'Invalid determinants: ' . join('.', @baddets) if @baddets;
            if (@dets == 1) {
                ++$singles{$dets[0]};
                ++$mults{$dets[0]};
            }
            elsif (@dets > 1) {
                foreach my $d (@dets) {
                    ++$mults{$d};
                }
                push @blends, join '.', @dets;
            }

            my $hasp = !!($fq eq '+');
            my $haso = !!($fq eq 'o');
            my $hasu = !!($fq eq 'u');
            my $hasm = !!($fq eq '-');
            my $hasn = !($hasp || $haso || $hasu || $hasm);

            $FQ{p} += $hasp;
            $FQ{o} += $haso;
            $FQ{u} += $hasu;
            $FQ{m} += $hasm;
            $FQ{n} += $hasn;

            # wrong, only if a single node matches M
            if (grep { /^M$/ } @dets) {
                $MQ{p} += $hasp;
                $MQ{o} += $haso;
                $MQ{u} += $hasu;
                $MQ{m} += $hasm;
                $MQ{n} += $hasn;
            }

            if ($hasW || $hasD) {
                $WDQ{p} += $hasp;
                $WDQ{o} += $haso;
                $WDQ{u} += $hasu;
                $WDQ{m} += $hasm;
                $WDQ{n} += $hasn;
            }

            $two =~ /^ (2?) $/x or die "Bad (2) column: $two";
            $twos += !!$1;

            my @conts = split /,/, $content;
            my @badconts = grep { !$contents{$_} } @conts;
            die 'Invalid contents: ' . join('.', @badconts) if @badconts;
            foreach my $cont (@conts) {
                ++$ccounts{$cont};
            }

            $pop =~ /^ (P?) $/x or die "Bad Pop column: $pop";
            $popular += !!$1;

            ++$Zf if $zscore;
            $ZSum += $zscore||0;

            my @specs = map { $specialsynonyms{$_} || $_ }
                        split /,/, $specscores;
            my @badspecs = grep { !exists($special{$_}) } @specs;
            die 'Invalid special scores: ' . join('.', @badspecs) if @badspecs;
            foreach my $spec (@specs) {
                ++$sc{$spec};
                $raw6 += $special{$spec} > 0;
                $weight6 += $special{$spec};
            }

            $sminus += $hasS && $hasm;
        }
    }

    my $ZEst = zest_table($Zf);

    my $W_D = $W + $D;

    my %qx = (p=>'+', m=>'-', n=>'none');
    my $fq = join "\n", map {
            sprintf("%4s    = %2d    = %2d    = %2d",
                    $qx{$_}||$_, $FQ{$_}, $MQ{$_}, $WDQ{$_});
        } qw(p o u m n);

    my $blends = join "\n", @blends;

    my $singles = join "\n", map {
            sprintf("%8s = %d", $_, $singles{$_}||0)
        } @determinants;

    my $ccounts = join "\n", map {
            sprintf("%8s = %d", $_, $ccounts{$_}||0)
        } @contents;

    my $approach = join "\n", map {
        sprintf("%5s   %s", $romans[$_], join('.', @{$approaches[$_]}));
    } 1..$#approaches;

    my $speciallist2 = join "\n", map {
        sprintf("%4s    = %2d x%1d  %2d x%1d", $_, $sc{$_}, $special1{$_},
                                               $sc{$_ . '2'}, $special2{$_})
    } qw(DV INC DR FAB);
    my $speciallist1 = join "\n", map {
        sprintf("%4s    = %2d x%1d", $_, $sc{$_}, $special1{$_})
    } qw(ALOG CON);

    my $lambdadenom = $respnum - $mults{F};
    my $lambda = $mults{F} / $lambdadenom;
    my $WSumC = 0.5 * $mults{FC} + $mults{CF} + 1.5 * $mults{C};
    my ($EBmin, $EBmax);
    if ($WSumC > $mults{M}) {
        $EBmin = $mults{M};
        $EBmax = $WSumC;
    }
    else {
        $EBmin = $WSumC;
        $EBmax = $mults{M};
    }
    my $EA = $mults{M} + $WSumC;
    my $EBPer = (
           $lambda < 1.0
        &&
           (   (4 <= $EA && $EA <= 10 && $EBmax >= $EBmin + 2)
            || ($EA > 10 && $EBmax >= $EBmin + 2.5)
           )
        && $EBmin != 0
        ) ? $EBmax / $EBmin : '-';
    #sum(m) and sumFM only include whole determinant components
    #sumC', sumT, sumV, sumY include sub-determinants
    my %detsums;
    my $sumshade = 0;
    foreach my $subdet ("C'", qw(T V Y)) {
        my $sum = 0;
        foreach my $det (grep { /$subdet/ } @determinants) {
            $sum += $mults{$det};
        }
        $detsums{$subdet} = $sum;
        $sumshade += $sum;
    }
    my $eb = $mults{FM} + $mults{m};
    my $es = $eb + $sumshade;
    my $adjes = $es - ($mults{m} > 1 ? $mults{m} - 1 : 0) 
                    - ($detsums{Y} > 1 ? $detsums{Y} - 1 : 0);

    my $d = d_table($EA - $es);
    my $adjd = d_table($EA - $adjes);

    #print "DEBUG dtable:\n", join "\n", map {
    #      "$_ to " . ($_ + 2) . ": "
    #    . d_table($_) . " and " . d_table($_ + 2)
    #} (13, 10.5, 8, 5.5, 3, -2.5, -5, -7.5, -10, -12.5, -15);

    my $abartay = 2*$sc{AB} + $ccounts{Art} + $ccounts{Ay};
    my $lv2 = $sc{DV2} + $sc{INC2} + $sc{DR2} + $sc{FAB2};

    my $CFpC = $mults{CF} + $mults{C};
    my $Afr = ($resps[8] + $resps[9] + $resps[10])
            / (  $resps[1] + $resps[2] + $resps[3] + $resps[4]
               + $resps[5] + $resps[6] + $resps[7]);
    my $numblends = scalar @blends;

    my $XAp = ($FQ{p} + $FQ{o} + $FQ{u}) / $respnum;
    my $WDAp = ($WDQ{p} + $WDQ{o} + $WDQ{u}) / ($W + $D);
    my $Xminusp = $FQ{m} / $respnum;
    my $Xplusp = $FQ{p} / $respnum;
    my $Xup =  $FQ{u} / $respnum;

    my $Zd = $ZEst eq '-' ? '-' : ($ZSum - $ZEst);

    my $human = $ccounts{H}  + $ccounts{'(H)'}
              + $ccounts{Hd} + $ccounts{'(Hd)'};
    my $isol_index = (  $ccounts{Bt} + 2*$ccounts{Cl} + $ccounts{Ge}
                      + $ccounts{Ls} + 2*$ccounts{Na}  ) / $respnum;

    my $frrf = $mults{Fr} + $mults{rF};
    my $egocentricity = (3*$frrf + $twos) / $respnum;
    my $anxy = $ccounts{An} + $ccounts{Xy};
    my $hhdhd = $ccounts{'(H)'} + $ccounts{Hd} + $ccounts{'(Hd)'};

    my $colshade = !! grep {   / C ([^']|$) .* (T|Y|V|C') /x
                            || / (T|Y|V|C') .* C ([^']|$) /x } @blends;
    my $pti = ($XAp < 0.70 && $WDAp < 0.75)
            + ($Xminusp > 0.29)
            + ($lv2 > 2 && $sc{FAB2} > 0)
            + (   ($respnum < 17 && $weight6 > 12)
               || ($respnum > 16 && $weight6 > 17) )
            + ($MQ{m} > 1 || $Xminusp > 0.40);
    my $depi = ($mults{FV} + $mults{VF} + $mults{V} > 0 || $mults{FD} > 2)
             + ($colshade || $S > 2)
             + ( ($egocentricity > 0.44 && $frrf == 0) || $egocentricity < 0.33 )
             + ($Afr < 0.46 || $numblends < 4)
             + ($sumshade > $eb || $detsums{"C'"} > 2)
             + ($sc{MOR} > 2 || $abartay > 3)
             + ($sc{COP} < 2 || $isol_index > 0.24);
    my $cdi = ($EA < 6 || $adjd < 0)
            + ($sc{COP} < 2 && $sc{AG} < 2)
            + ($WSumC < 2.5 || $Afr < 0.46)
            + ($ap{p} > $ap{a}+1 || $ccounts{H} < 2)
            + ($detsums{T} > 1 || $isol_index > 0.24 || $ccounts{Fd} > 0);
    my $scons = ($mults{FV} + $mults{VF} + $mults{V} + $mults{FD} > 2)
              + $colshade
              + ($egocentricity > 0.44 || $egocentricity < 0.31)
              + ($sc{MOR} > 3)
              + ($Zd > 3.5 || $Zd < -3.5)
              + ($es > $EA)
              + ($CFpC > $mults{FC})
              + ($Xplusp < 0.70)
              + ($S > 3)
              + ($popular < 3 || $popular > 8)
              + ($ccounts{H} < 2)
              + ($respnum < 17);
    my $hdad = $ccounts{Hd} + $ccounts{Ad};
    my $haratio = ($hdad != 0) ? ($ccounts{H} + $ccounts{A}) / $hdad : 0;
    my $hvisubs = ($Zf > 12)
                + ($Zd > 3.5)
                + ($S > 3)
                + ($ccounts{H} + $hhdhd > 6)
                + (  $ccounts{'(H)'} + $ccounts{'(A)'} + $ccounts{'(Hd)'}
                   + $ccounts{'(Ad)'} > 3 )
                + ($haratio < 4)
                + ($ccounts{Cg} > 3);
    my $hvi = ($detsums{T} == 0 && $hvisubs >= 4) ? 'Yes' : 'No';
    my $obssum1_4 = ($Dd > 3)
                  + ($Zf > 12)
                  + ($Zd > 3.0)
                  + ($popular > 7);
    my $obssum1_5 = $obssum1_4 + ($FQ{p} > 1);
    my $obs = (  $obssum1_5 == 5
              || ($obssum1_4 >= 2 && $FQ{p} > 3)
              || ($obssum1_5 >= 3 && $Xplusp > 0.89)
              || ($FQ{p} > 3 && $Xplusp > 0.89) ) ? 'Yes' : 'No';

    foreach (@special0) {
        $sc{$_} = sprintf("%2d", $sc{$_}||0);
    }

    return <<END;
Location Features:
Zf      = $Zf
ZSum    = $ZSum
ZEst    = $ZEst

W       = $W
D       = $D
W+D     = $W_D
Dd      = $Dd
S       = $S

DQ
+       = $DQ{p}
o       = $DQ{o}
v/+     = $DQ{vp}
v       = $DQ{v}

Form Quality:
        FQx     MQual   W+D
$fq

Determinants:
Blends
$blends

Single
$singles

(2)     = $twos

Contents:
$ccounts

Approach:
$approach

Special Scores:
          Lv1   Lv2
$speciallist2
$speciallist1
Raw Sum6  = $raw6
Wgtd Sum6 = $weight6

AB      = $sc{AB}     GHR = $sc{GHR}
AG      = $sc{AG}     PHR = $sc{PHR}
COP     = $sc{COP}     MOR = $sc{MOR}
CP      = $sc{CP}     PER = $sc{PER}
                 PSV = $sc{PSV}

Ratios, Percentages, and Derivations
R = $respnum          L      = $mults{F} / $lambdadenom = $lambda
EB = $mults{M} : $WSumC    EA     = $EA    EBPer = $EBPer
eb = $eb : $sumshade    es     = $es      D = $d
                Adj es = $adjes  Adj D = $adjd

FM = $mults{FM}    SumC' = $detsums{"C'"}   SumT = $detsums{T}
m  = $mults{m}    SumV  = $detsums{V}   SumY = $detsums{Y}

a:p     = $ap{a} : $ap{p}
Ma:Mp   = $Map{a} : $Map{p}
2AB+Art+Ay = $abartay
MOR     = $sc{MOR}

Sum6    = $raw6
Lv2     = $lv2
WSum6   = $weight6
M-      = $MQ{m}
Mnone   = $MQ{n}


FC:CF+C     = $mults{FC} : $CFpC
Pure C      = $mults{C}
SumC':WSumC = $detsums{"C'"} : $WSumC
Afr         = $Afr
S           = $S
Blends:R    = $numblends : $respnum
CP          = $sc{CP}

XA%     = $XAp
WDA%    = $WDAp
X-%     = $Xminusp
S-      = $sminus
P       = $popular
X+%     = $Xplusp
Xu%     = $Xup

Zf      = $Zf
W:D:Dd  = $W : $D : $Dd
W:M     = $W : $mults{M}
Zd      = $Zd
PSV     = $sc{PSV}
DQ+     = $DQ{p}
DQv     = $DQ{v}

COP = $sc{COP} AG = $sc{AG}
GHR:PHR     = $sc{GHR} : $sc{PHR}
a:p         = $ap{a} : $ap{p}
Food        = $ccounts{Fd}
SumT        = $detsums{T}
Human Cont  = $human
Pure H      = $ccounts{H}
PER         = $sc{PER}
Isol Indx   = $isol_index


3r+(2)/R        = $egocentricity
Fr+rF           = $frrf
SumV            = $detsums{V}
FD              = $mults{FD}
An+Xy           = $anxy
MOR             = $sc{MOR}
H:(H)+Hd+(Hd)   = $ccounts{H} : $hhdhd

PTI     = $pti
DEPI    = $depi
CDI     = $cdi
S-Con   = $scons
HVI     = $hvi
OBS     = $obs
END
}

sub d_table {
    my ($input) = @_;
    die "EA - es is too small for D table" if $input < -15.0;
    die "EA - es is too big for D table" if $input > 15.0;
    if ($input >= 0) {
        return int(($input - 0.25) / 2.5)
    }
    else {
        return -int((-$input - 0.25) / 2.5);
    }
}

sub zest_table {
    my ($zf) = @_;
    die "Zf $zf is too low" if $zf < 1;
    die "Zf $zf is too high" if $zf > 50;
    my @zest_table = (
         '-',    2.5,   6.0,  10.0,  13.5,  17.0,  20.5,  24.0,  27.5,  31.0,
         34.5,  38.0,  41.5,  45.5,  49.0,  52.5,  56.0,  59.5,  63.0,  66.5,
         70.0,  73.5,  77.0,  81.0,  84.5,  88.0,  91.5,  95.0,  98.5, 102.5,
        105.5, 109.5, 112.5, 116.5, 120.0, 123.5, 127.0, 130.5, 134.0, 137.5,
        141.0, 144.5, 148.0, 152.0, 155.5, 159.0, 162.5, 166.0, 169.5, 173.0
    );
    return $zest_table[$zf - 1] or die "Zf $zf out of bounds";
}
