use strict;
use warnings;

my ($in_tables_dir, $out_tables_file_full_path) = @ARGV;

if($in_tables_dir =~ m/(.*)\/$/)
{
	# remove redundant '/' if exists:
   	$in_tables_dir = $1;
}

my @all_tables_files = glob($in_tables_dir . '/chop_tables*');
open (my $out, ">", $out_tables_file_full_path) or die "could not open $out_tables_file_full_path for writing";

# ~/long_indel/chop_tables/chop_tables_r_param_0.3_basic_mu_0.04_basic_gamma_0.99_max_indel_length_50_length_of_anc_1000000_t_0.25.txt
# ~/long_indel/chop_tables_INDELible/chop_tables_A_param_1.1_IR_0.01_max_indel_length_50_length_of_anc_1000000_t_0.25.txt
foreach my $chop_tables_file (@all_tables_files)
{
	my $curr_r_or_A_param;
	my $curr_basic_mu_or_IR;
	my $curr_basic_gamma = 0.0;
	my $curr_max_indel_length;
	my $curr_t;
	if ($chop_tables_file =~ m/(r_param|A_param)_(\d(\.\d+)?)_(basic_mu|IR)_(\d(\.\d+)?)_(basic_gamma_(\d(\.\d+)?)_)?max_indel_length_(\d+).*?_t_(\d(\.\d+)?)/)
	{
		$curr_r_or_A_param = $2;
		$curr_basic_mu_or_IR = $5;
		if (defined $8)
		{
			$curr_basic_gamma = $8;
		}
		$curr_max_indel_length = $10;
		$curr_t = $11;
	}
	
	print $out "$chop_tables_file\t$curr_r_or_A_param\t$curr_basic_mu_or_IR\t$curr_basic_gamma\t$curr_max_indel_length\t$curr_t\n";
}
close ($out);