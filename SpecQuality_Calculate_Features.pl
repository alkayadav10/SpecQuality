#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum max min);
use constant H=>1.007825;
use XML::Twig;
use MIME::Base64;
use Compress::Zlib;

# Program to read MS/MS spectra files (MGF, mzML)
# and process each spectrum to assign a quality score

# INPUT: MGF & mzML
# 21 Aug 2025: added mzML
# # 21:41 14 October 2024
# Updated to version 6 Date : 29 April 2025
# OUTPUT : Scan,SQS (Spectral Quality Score)


######### INPUT PARAMETERS ####################
my $input_file=shift;
my $msgf =shift;

my $tol = 0.5;
my $tolunit = 'Da'; # Da or ppm
# Define a noise threshold for S/N calculation
my $noise_threshold = 5;  # Threshold in terms of % normalized intensity
my $peakcutoff= 7;
###############################################

my $t0=time();
chomp $input_file;

print "Starting SpecQuality analysis for: $input_file\n";

# Define amino acid masses (in Da)
my %amino_acid_masses = (
    'A' => 71.03711,   # Alanine
    'C' => 103.00919,  # Cysteine
    'D' => 115.02694,  # Aspartic acid
    'E' => 129.04259,  # Glutamic acid
    'F' => 147.06841,  # Phenylalanine
    'G' => 57.02146,   # Glycine
    'H' => 137.05891,  # Histidine
    'I' => 113.08406,  # Isoleucine
    'K' => 128.09496,  # Lysine
    'L' => 113.08406,  # Leucine
    'M' => 131.04049,  # Methionine
    'N' => 114.04293,  # Asparagine
    'P' => 97.05276,   # Proline
    'Q' => 128.05858,  # Glutamine
    'R' => 156.10111,  # Arginine
    'S' => 87.03203,   # Serine
    'T' => 101.04768,  # Threonine
    'V' => 99.06841,   # Valine
    'W' => 186.07931,  # Tryptophan
    'Y' => 163.06333,  # Tyrosine
);

my $out=$input_file;
$out=~s/\.(mgf|mzML)$/_$1_features.tsv/i;

my $PSM_SQS=$input_file;
$PSM_SQS=~s/\.(mgf|mzML)$/_$1_PSM_features.tsv/i;
open (OUT,">$out") or die $!;
print OUT "scan\tmz\tz\tRT\tPeakCount\tSignalPeaksCount\tSNR\tIntensity\tPeakDensity\tGoodDiffFraction\tComplementsFraction\tIsotopePeaks\tNeutralLossPeaks\tAverageRelativeIntensity\tSQS_GM\n";

my %sp = parse_ms_file($input_file, $peakcutoff, $tol, $tolunit, $noise_threshold);
my $t1=time();
print "-> File parsing completed in : ", $t1-$t0," seconds (",($t1-$t0)/60," minutes)\n";
print "   Total spectra loaded: ", scalar(keys %sp), "\n";


#normalize each factor
print "Normalizing 'complements_fraction' feature...\n";
my $sref = fix_complements_fraction(\%sp);
%sp =%$sref;

# my %spnew;
# foreach my $scan_id (keys %sp) {
#     # This copies the inner hash for each scan_id
#     $spnew{$scan_id} = { %{ $sp{$scan_id} } };
# }

# print "Performing range normalization on all features\n";
# %spnew = range_normalize(\%spnew,\%minmax);
print "Calculating Spectral Quality Score (SQS) for each spectrum...\n";
%sp = calculate_SQS(%sp);

print "-> SQS calculation completed for ", scalar keys %sp," spectra\n";

#PRINT in OUTPUT FILE
print "Writing SQS results to output file: $out\n";
foreach my $scan(keys %sp)
{
	my $features = $sp{$scan};
    my ($mz,$z,$RT) = @{$sp{$scan}{COL}};
    my $original_title = $sp{$scan}{original_title} // $scan;  # Use original title if available, fallback to the key
    
    # print "Before:\t",$sp{$scan}{peak_count},"\n";
    # print "After:\t",$spnew{$scan}{peak_count},"\n";<>;

    # Output results with new features, using original title
    print OUT join("\t", 
            $original_title,  # scan (original title)
            $mz,              # mz
            $z,               # z
            $RT,              # RT
            
            # RAW detailed features
            $features->{peak_count},
            $features->{signal_peaks},
            $features->{snr},
            $features->{intensity},
            $features->{peak_density},
            $features->{good_diff_fraction},
            $features->{complements_fraction},
            $features->{isotope_peaks},
            $features->{neutral_loss_peaks},
            $features->{ari},
            $features->{sqs},
        ). "\n";	
}
close OUT;
my $t2=time();
print "SQS calculated and Features printed in : ", $t2-$t1," seconds (",($t2-$t1)/60," minutes)\n";

if((scalar @ARGV>=1)&&(defined $msgf))
{
print "Merging SQS results with PSM file: $msgf\n";


open MSGF,$msgf or die $!;
open MSGFOUT,">$PSM_SQS" or die $!;
my $f=0; #header flag

#Map scans and print SQS results merged with PSM file
my $c=0; #non matched
my $p=0;#matched spec
my $t=0;#total spec

while(<MSGF>)
{
	chomp $_;
	if ($f==0)
	{
		$f=1;
		#print header #"scan\tmz\tz\tRT\tpeaks\tSignalPeaks\tSNRatio\tMS2Int\tnormMS2Int\tMS2SumRelInt\tMS1Int\tnormMS1Int\tSQS\n";
		print MSGFOUT "$_\tPeakCount\tSignalPeaksCount\tSNR\tIntensity\tPeakDensity\tGoodDiffFraction\tComplementsFraction\tIsotopePeaks\tNeutralLossPeaks\tAverageRelativeIntensity\tSQS\tLabel\n";
        next;
	}
    $t++;
	my @arr=split/\t/,$_;
	my $scan_title=$arr[3];
    my $charge = $arr[8]; # Assuming charge is in the 9th column
    my $qval=$arr[16];

    # Find the correct unique key in %sp that matches the scan title and charge
    my $matched_key;
    foreach my $unique_key (keys %sp) {
        # Reconstruct a partial key from PSM data to find the match
        if ($unique_key =~ /^\Q$scan_title\E.*_z_\Q$charge\E$/) {
            $matched_key = $unique_key;
            last;
        }
    }

	if (defined $matched_key && exists $sp{$matched_key})
	{
        $p++;
        my $features = $sp{$matched_key};
        # my $final = $spnew{$scan};
        #my @SQS=@{$sp{$scan}{SQS}};
        my $label ='NA';
        if ($qval<=0.01)
        {
            $label = 'Excellent';
        }
        elsif ($qval <= 0.05)
        {
            $label = 'Good';
        }
        elsif ($qval <= 0.1)
        {
            $label = 'Average';
        }
        elsif ($qval > 0.1)
        {
            $label = 'Poor';
        }

        print MSGFOUT "$_\t",join("\t",
         # New detailed features
            $features->{peak_count},
            $features->{signal_peaks},
            $features->{snr},
            $features->{intensity},
            $features->{peak_density},
            $features->{good_diff_fraction},
            $features->{complements_fraction},
            $features->{isotope_peaks},
            $features->{neutral_loss_peaks},
            $features->{ari},   
            $features->{sqs},
            $label
        )."\n";
        delete $sp{$matched_key};#delete scan to avoid repeats	
    }
    else
    {
        $c++;
    }
}

print "MSGF+ scans(Total):$t, Matched: $p and Unmatched: $c scans\n";

close MSGF;
close MSGFOUT;
}

my $t3=time();
print "Total time taken : ", $t3-$t0," seconds (",($t3-$t0)/60," minutes)\n";

exit;

sub calculate_spectral_quality_features {
    my ($peaks_ref, $parent_MH, $intensity, $z, $tol, $tolunit,$noise_threshold) = @_;

    my @peaks = @$peaks_ref;
    # my $ticlog = log10($tic); 
    $intensity = log10($intensity);
    
    # Peak Count (Npeaks)
    my $peak_count = scalar @peaks;

    # First, find the maximum intensity for spectrum normalization
    my $max_intensity = 0;
    foreach my $peak (@peaks)
	{
        my ($mz, $int) = @$peak;
        if ($int > $max_intensity)
		{
            $max_intensity = $int;
        }
    }

    # Signal Peaks (SP)
    my $signal_peaks = 0;
    
    my $total_signal_intensity = 0; #sum(map { $_->[1] } grep { $_->[1] >= $max_intensity * $noise_threshold } @peaks);
    my $total_noise_intensity = 1; #sum(map { $_->[1] } grep { $_->[1] < $max_intensity * $noise_threshold } @peaks);
    
    my $minmz = $parent_MH;
    my $maxmz = 0;
    # Now calculate normalized intensities and assess quality for SNR
    foreach my $peak (@peaks)
	{
        my ($mz, $int) = @$peak;
        my $normalized_intensity = ($int*100)/ $max_intensity;
        
        # Signal-to-noise ratio calculation using normalized intensity
        if ($normalized_intensity >= $noise_threshold)
		{
            $signal_peaks++;
            $total_signal_intensity += $int;
        }
        else
        {
            $total_noise_intensity += $int;
        }

        if ($minmz > $mz)
		{
            $minmz = $mz;
        }
        if ($maxmz< $mz)
        {
            $maxmz = $mz;
        }
    }
    
    # Signal-to-Noise Ratio (SNR)
    my $snr = log10(($total_signal_intensity)/($total_noise_intensity));
    if ($snr <=0)
    {
        $snr = 0.01;
    }
    #print "$title\t$snr\n" if $snr<0;
    #print "SNR=$snr\n";<STDIN>;
        
    # Peak Density
    my $mz_range = $maxmz - $minmz;
    my $peak_density = $peak_count*100 / $mz_range;
    
    # Good-Diff Fraction and Complements
    my $good_diffs = 0;
    my $total_diffs = 0;
    my $complements = 0;
	for (my $i=0;$i<@$peaks_ref - 1;$i++)
	{
		my ($mz_i, $intensity_i) = @{$peaks_ref->[$i]};
		#print"mz=$mz_i\tint=$intensity_i\n";
		for (my $j= ($i + 1);$j<@$peaks_ref;$j++)
		{			
			my ($mz_j, $intensity_j) = @{$peaks_ref->[$j]};
			my $diff = abs($mz_i - $mz_j);
			#print"\tmz=$mz_j\tint=$intensity_j\tmzdiff=$diff\n";
			my $good_diff_found = 0;
			# Total Diff within specified range
			if ($diff >= 56 && $diff <= 187) # Within AA mass range
			{  
				$total_diffs += $intensity_i + $intensity_j ; #denominator for normalizing
				
				#print"\t\ttotalDiffscumul=$total_diffs\n";
				# Good-Diff Fraction
				
				foreach my $mass (values %amino_acid_masses)
				{
					if (($tolunit eq 'Da') && (abs($diff - $mass) <= $tol))
					{
						$good_diff_found = 1;
						#print"\t\t\tgooddiff=$mass\n";
						last;  # Exit loop as soon as a match is found
					}
                    elsif ($tolunit eq 'ppm')
                    {
                        my $ppmtol = ppm_to_da($tol,$diff);
                        if ((abs($diff - $mass) <= $ppmtol))
                        {
                            $good_diff_found = 1;
                            last;  # Exit loop as soon as a match is found
                        }
                    }
				}
			}

			if ($good_diff_found)
			{
				$good_diffs += $intensity_i + $intensity_j;
				#print"\t\tgoodDiffscumul=$good_diffs\n";
			}

			# Complements
			my $sum = $mz_i + $mz_j;
			my $diff_from_parent = abs($parent_MH - $sum);
			if (($tolunit eq 'Da') && ($diff_from_parent <= $tol))
			{
				$complements += $intensity_i + $intensity_j;
				#print"\t\tcomplements-cumul=$complements\n";
			}
            elsif ($tolunit eq 'ppm')
			{
                my $ppmtol = ppm_to_da($tol,$diff_from_parent);
                if (abs($diff_from_parent) <= $ppmtol)
				{
                    $complements += $intensity_i + $intensity_j;
				    #print"\t\tcomplements-cumul=$complements\n";
                }
			}
			#<>;
		}
		#print"\n";
	}

    # Calculate Good-Diff Fraction, avoid zeroes
    my $good_diff_fraction = $total_diffs > 0 ? ( (100 * $good_diffs) / $total_diffs) : 0.01;
    if ($good_diff_fraction <= 0)
    {
        $good_diff_fraction = 0.01;
    }
	
	# Calculate Complements Fraction, avoid zeroes
    if ($complements <= 0)
    {
        $complements = 0.01;
    }
    my $complements_fraction = log10($complements/$intensity);
    #my $complements_fraction = $complement_count / (($peak_count * ($peak_count - 1) / 2) || 1);
    
    # # 8. Intensity Balance
    # my @intensity_bands = (0) x 10;
    # my $min_mz = $peaks[0][0];
    # #print "MIN=",$min_mz."\n";<>;
    # my $max_mz = max(map { $_->[0] } @peaks);
    # my $band_width = ($max_mz - $min_mz) / 10;
    
    # for my $peak (@peaks) {
    #     my $band_index = int(($peak->[0] - $min_mz) / $band_width);
    #     $band_index = 9 if $band_index > 9;
    #     $intensity_bands[$band_index] += $peak->[1];
    # }
    
    # my $lowest_7_bands_intensity = sum(@intensity_bands[0..6]);
    # my $highest_2_bands_intensity = sum(@intensity_bands[8..9]);
    # my $intensity_balance = $highest_2_bands_intensity - $lowest_7_bands_intensity;
    # #print "mz=$precursor_mz\tint=$intensity ($tic)\tRaw intbal=$intensity_balance\t\n";#if $intensity_balance>$tic;
    # $intensity_balance = log10($tic-$intensity_balance);
    # $intensity_balance = log10($tic-$intensity_balance)/$tic;
    #print "\t TIC norm intbal=$intensity_balance\t\n";

    #$intensity_balance = 1 - abs(log($intensity_balance)) / 2;
    #print "\t log norm intbal=$intensity_balance\t\n";<>;
#    print "\tmz=$precursor_mz\tint=$intensity ($tic)\tintbal=$intensity_balance\tIB=$ib\n";<>;
    
    # Isotope and Neutral Loss Peaks
    my $isotope_peaks = 0;
    my $neutral_loss_peaks = 0;
    
    for my $i (0..$#peaks-1)
    {
        # Isotope peaks
        if (($tolunit eq 'Da') && (abs($peaks[$i+1][0] - $peaks[$i][0] - 1.0) < $tol)) {
            $isotope_peaks++;
        }
        elsif($tolunit eq 'ppm')
        {
            my $ppmtol = ppm_to_da($tol,$peaks[$i][0]);
            if (abs($peaks[$i+1][0] - $peaks[$i][0] - 1.0) < $ppmtol)
            {
                $isotope_peaks++;
            }
        }
        
        # Neutral loss peaks
        my @neutral_loss_masses = (17, 18, 28);
        for my $loss (@neutral_loss_masses)
        {
            if (($tolunit eq 'Da') && (abs($peaks[$i+1][0] - $peaks[$i][0] - $loss) < $tol))
            {
                $neutral_loss_peaks++;
            }
            elsif ($tolunit eq 'ppm')
            {
                my $ppmtol = ppm_to_da($tol, $loss);
                if (abs($peaks[$i+1][0] - $peaks[$i][0] - $loss) < $ppmtol)
                {
                    $neutral_loss_peaks++;
                }
            }
        }
    }
    
    # Avoid numerator to be zero
    if ($isotope_peaks > 0)
    {
        $isotope_peaks = 100*$isotope_peaks/$peak_count;
    }
    else
    {
        $isotope_peaks = 0.01;
    }

    if ($neutral_loss_peaks > 0)
    {
        $neutral_loss_peaks = 100*$neutral_loss_peaks/$peak_count;
    }
    else
    {
        $neutral_loss_peaks = 0.01;
    }

    # Number of Dominant Peaks (NDP)
    #my @sorted_peaks = sort { $b->[1] <=> $a->[1] } @peaks;
    #my $total_intensity = $tic;
    # my $cumulative_intensity = 0;
    # my $ndp = 0;
    
    # #for my $peak (@sorted_peaks) {
    # for my $peak (@peaks) 
    # {
    #     $ndp++ if $peak->[1] >= 0.5 * $max_intensity;
    # }
    # $ndp = $ndp *100 / $precursor_mz;
    
    # Average of Relative Peak Intensities (ARI)
    my $ari = $intensity*$signal_peaks/$peak_count;

    # print "PREC=$precursor_mz\tINT=$intensity\n";
    # print "INSIDE\nPeak count: $peak_count\n";
    # print "Signal peaks: $signal_peaks\n";
    # print "SNR: $snr\n";
    # print "MS1 intensity: ",$intensity,"\n";    
    # print "TIC: ",$tic,"\t",$tic,"\n";
    # print "Peak density: $peak_density\n";
    # print "Good Diff Fraction: $good_diff_fraction\n";
    # print "Complement fraction: $complements_fraction\n";
    # print "Intensity balance: $intensity_balance\n";
    # print "Isotope peaks: $isotope_peaks\n";
    # print "Neutral loss peaks: $neutral_loss_peaks\n";
    # print "Iso+Neutral loss peaks: ", $isotope_peaks+$neutral_loss_peaks,"\n";
    # print "Number of dominant peaks (NDP): $ndp\n";
    # print "Average of relative peak intensities (ARI): ",$ari,"\n";<STDIN>;
    
    return {
        peak_count => $peak_count,
        signal_peaks => $signal_peaks,
        snr => $snr,
        intensity => $intensity,
        peak_density => $peak_density,
        good_diff_fraction => $good_diff_fraction,
        complements_fraction => $complements_fraction,
        isotope_peaks => $isotope_peaks,
        neutral_loss_peaks => $neutral_loss_peaks,
        # ndp => $ndp,
        ari => $ari,
        sqs => undef,
    };
}

# Calculate quality features for a given spectrum
sub calculate_SQS {
    my (%sp) = @_;
    
    # Features to include in SQS calculation
    my @features_to_use = qw(
        peak_count 
        signal_peaks 
        snr
        intensity 
        peak_density 
        good_diff_fraction 
        complements_fraction 
        isotope_peaks 
        neutral_loss_peaks 
        ari
    );
    
    
    # my $sum=0;
    foreach my $scan(keys %sp)
    {
        # Calculate geometric mean
        my $product = 1;
        $product *= $sp{$scan}{$_} for @features_to_use;
        my $geometric_mean = $product ** (1 / scalar(@features_to_use));
        
        $sp{$scan}{sqs} = $geometric_mean;
        
    }
    return %sp;
}


# # Subroutine to convert ppm tolerance to Da
sub ppm_to_da {
    my ($ppm, $mass) = @_;
    # print "\tppm=$ppm\tmass=$mass\n";
    my $ppmtol = ($ppm * $mass) / 1_000_000;
    # print "\tppmtol=$ppmtol\n";
    return ($ppmtol);
}

# Function for log10 calculator 
sub log10
{ 
    my ($n) = @_;      
    # using pre-defined log function 
    my $log10 = log($n) / log(10);
    return ($log10);
} 

sub fix_complements_fraction{
    my ($spref) = @_; 
    my %sp=%$spref;
    #print "\tFixing complements fraction\n";
    my $min=100000000;
    my $max =0;
    
    foreach my $spectrum (keys %sp)
    {
        $min = $sp{$spectrum}{complements_fraction} if $min > $sp{$spectrum}{complements_fraction};
        $max = $sp{$spectrum}{complements_fraction} if $max < $sp{$spectrum}{complements_fraction};
    }

    if($min< 0)
    {
        foreach my $spectrum (keys %sp)
        {
            $sp{$spectrum}{complements_fraction} += -($min)+0.01;
        }
    }
    else
    {
        return (\%sp);  # Return the original spectrum hash
    }
    # print "\tFixed complements fraction\n";
    # print "\tMin=$minmax{complements_fraction}[0]\tMax=$minmax{complements_fraction}[1]\n";

    return (\%sp);
}

sub range_normalize {
    my ($spref, $ref) = @_;
    my %spnew=%$spref;
    my %minmax=%$ref;

    foreach my $spectrum(keys %spnew)
    {
        my @features = qw(
        peak_count 
        signal_peaks 
        snr
        intensity 
        peak_density 
        good_diff_fraction 
        complements_fraction 
        isotope_peaks 
        neutral_loss_peaks 
        ari
        );

        foreach my $feat(@features)
        {
            # print "Feature=$feat\tMin=$minmax{$feat}[0]\tMax=$minmax{$feat}[1]\n";
            $spnew{$spectrum}{$feat} = 100*($spnew{$spectrum}{$feat})/($minmax{$feat}[1] - $minmax{$feat}[0]);
            if ($spnew{$spectrum}{$feat} == 0)
            {
                #print "Scan=$spectrum\tFeat=$feat\tValue=$sp{$spectrum}{$feat}\n";
                $spnew{$spectrum}{$feat} = 0.01;
            }
        }
    }
    return %spnew;
}

sub parse_ms_file {
    my ($file_path, $peakcutoff, $tol, $tolunit, $noise_threshold) = @_;

    if ($file_path =~ /\.mgf$/i) {
        return parse_mgf($file_path, $peakcutoff, $tol, $tolunit, $noise_threshold);
    } elsif ($file_path =~ /\.mzML$/i) {
        return parse_mzml($file_path, $peakcutoff, $tol, $tolunit, $noise_threshold);
    } else {
        die "Unsupported file format: $file_path. Only .mgf and .mzML are supported.";
    }
}

use Fcntl 'O_RDONLY';

sub parse_mgf {
    my ($mgf_file, $peakcutoff, $tol, $tolunit, $noise_threshold) = @_;
    my %sp;
    
    # Use sysopen for faster file reading
    sysopen(my $MGF, $mgf_file, O_RDONLY) or die "Cannot open $mgf_file: $!\n";
    
    # Read entire file into memory for faster processing
    my $content = do { local $/; <$MGF> };
    close $MGF;
    
    my $count = 0;
    my $pos = 0;
    my $content_length = length($content);
    
    # Process the file in chunks to find spectra
    while ($pos < $content_length) {
        # Find the next spectrum
        my $start = index($content, "BEGIN IONS\n", $pos);
        last if $start == -1;
        
        my $end = index($content, "END IONS\n", $start);
        last if $end == -1;
        
        my $spectrum = substr($content, $start + 11, $end - $start - 11);
        $pos = $end + 8;  # Move past END IONS
        
        # Parse spectrum
        my ($title, $mz, $int, $rt, $z) = ("unknown_$count", 0, 1, 0, 0);
        my @current_peaks = ();
        my $in_peaks = 0;
        
        # Process spectrum lines
        foreach my $line (split(/\n/, $spectrum)) {
            next if $line =~ /^\s*$/;  # Skip empty lines
            
            if ($line =~ /^TITLE=(.*)/) {
                $title = $1;
            }
            elsif ($line =~ /^PEPMASS=([\d\.]+)(?:\s+([\d\.]+))?/) {
                $mz = $1;
                $int = $2 if defined $2;
            }
            elsif ($line =~ /^RTINSECONDS=([\d\.]+)/) {
                $rt = $1;
            }
            elsif ($line =~ /^CHARGE=(\d+)\+?/) {
                $z = $1;
            }
            elsif ($line =~ /^(\d+\.?\d*)\s+(\d+\.?\d*)/) {
                push @current_peaks, [$1, $2];
            }
        }
        
        # Process the spectrum if it has enough peaks
        if (@current_peaks >= $peakcutoff) {
            $count++;
            if ($count % 1000 == 0) {
                print "  ... processed $count spectra\n";
            }
            
            my $precursorMH = $mz * $z;
            my $features = calculate_spectral_quality_features(
                \@current_peaks, 
                $precursorMH, 
                $int, 
                $z, 
                $tol, 
                $tolunit, 
                $noise_threshold
            );

            my $unique_scan_key = sprintf("%s_mz_%.4f_z_%d", $title, $mz, $z);
            $sp{$unique_scan_key} = $features;
            $sp{$unique_scan_key}{COL} = [$mz, $z, $rt];
            $sp{$unique_scan_key}{original_title} = $title;
        }
    }
    
    print "  Successfully processed $count spectra\n";
    return %sp;
}

sub parse_mzml {
    my ($mzml_file, $peakcutoff, $tol, $tolunit, $noise_threshold) = @_;
    my %sp;
    my $count = 0;

    my $twig = XML::Twig->new(
        twig_handlers => {
            'spectrum' => sub { _process_spectrum(@_, \%sp, \$count, $peakcutoff, $tol, $tolunit, $noise_threshold) },
        }
    );
    $twig->parsefile($mzml_file);
    return %sp;
}

sub _process_spectrum {
    my ($twig, $spectrum_element, $sp_ref, $count_ref, $peakcutoff, $tol, $tolunit, $noise_threshold) = @_;

    my %spec_data;
    $spec_data{title} = $spectrum_element->att('id') || 'Unknown';

    my $ms_level_cv = $spectrum_element->first_child('cvParam[@accession="MS:1000511"]');
    return unless ($ms_level_cv && $ms_level_cv->att('value') == 2);

    my $rt_cv = $spectrum_element->first_child('scanList')->first_child('scan')->first_child('cvParam[@accession="MS:1000016"]');
    my $rt = $rt_cv ? $rt_cv->att('value') : 0;

    my $mz = 0;
    my $z = 0;
    my $intensity = 1;
    my $precursor = $spectrum_element->first_child('precursorList')->first_child('precursor');
    if ($precursor) {
        my $selected_ion = $precursor->first_child('selectedIonList')->first_child('selectedIon');
        if ($selected_ion) {
            my $mz_cv = $selected_ion->first_child('cvParam[@accession="MS:1000744"]');
            $mz = $mz_cv ? $mz_cv->att('value') : 0;

            my $charge_cv = $selected_ion->first_child('cvParam[@accession="MS:1000041"]');
            $z = $charge_cv ? $charge_cv->att('value') : 0;

            my $int_cv = $selected_ion->first_child('cvParam[@accession="MS:1000042"]');
            $intensity = $int_cv ? $int_cv->att('value') : 1;
        }
    }

    my $binary_data_arrays = $spectrum_element->first_child('binaryDataArrayList');
    return unless $binary_data_arrays;

    my ($mz_array, $int_array);
    foreach my $binary_data_array ($binary_data_arrays->children('binaryDataArray')) {
        if ($binary_data_array->first_child('cvParam[@accession="MS:1000514"]')) {
            $mz_array = $binary_data_array;
        } elsif ($binary_data_array->first_child('cvParam[@accession="MS:1000515"]')) {
            $int_array = $binary_data_array;
        }
    }

    return unless ($mz_array && $int_array);

    my $decode_binary = sub {
        my ($array) = @_;
        my $is_zlib = $array->first_child('cvParam[@accession="MS:1000574"]');
        my $is_64bit = $array->first_child('cvParam[@accession="MS:1000523"]');
        my $is_32bit = $array->first_child('cvParam[@accession="MS:1000521"]');

        my $binary_element = $array->first_child('binary');
        return [] unless ($binary_element && $binary_element->text);

        my $binary_text = $binary_element->text;
        my $decoded = decode_base64($binary_text);

        $decoded = uncompress($decoded) if $is_zlib;

        my $unpack_format = 'd*';
        if ($is_32bit) {
            $unpack_format = 'f*';
        }
        return [unpack($unpack_format, $decoded)];
    };

    my $mz_values = $decode_binary->($mz_array);
    my $int_values = $decode_binary->($int_array);

    my @peaks;
    for (my $i = 0; $i < @$mz_values; $i++) {
        push @peaks, [$mz_values->[$i], $int_values->[$i]];
    }

    if (scalar @peaks < $peakcutoff) {
        return;
    }

    $$count_ref++;
    if ($$count_ref % 1000 == 0) {
        print "  ... processed $$count_ref spectra\n";
    }

    my $precursorMH = ($mz * $z);
    my $features = calculate_spectral_quality_features(\@peaks, $precursorMH, $intensity, $z, $tol, $tolunit, $noise_threshold);

    # Create a more unique key to avoid overwriting spectra with same m/z but different charge
    my $unique_scan_key = sprintf("%s_mz_%.4f_z_%d", $spec_data{title}, $mz, $z);

    $sp_ref->{$unique_scan_key} = $features;
    $sp_ref->{$unique_scan_key}{COL} = [$mz, $z, $rt];
    $sp_ref->{$unique_scan_key}{original_title} = $spec_data{title};  # Store the original title

    $spectrum_element->purge();
}