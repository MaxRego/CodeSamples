#! /usr/bin/perl
#
######################################################################################################
#
#	Program: 	query_dmel.pl
#	Version:	1.0
#
#	Author:		Max Rego
#	Email:		mr255509@ohio.edu
#	Project:	Drosophila Melanogaster Orthologous Genes Database Query Script
#	
#	Description:	This program access the database built by DrosphilaMelanogaster.pl and prints
#			useful information to a datafile
#
#	Date:		3-14-2013
#
######################################################################################################

use DBI;
use DBD::mysql;
use Bio::SeqIO;
use Cwd;

#Connecting to database
my $ds = "DBI:mysql:dmel:localhost";
my $user = "";
my $passwd = "";
my $dbh = DBI -> connect($ds, $user, $passwd) || die "Cant Connect!\n";

#values from XML
my $jobID = $ARGV[0];
my $inFile = $ARGV[1];
my $upLen = $ARGV[2];
my $downLen = $ARGV[3];
my $speciesString = $ARGV[4];
my $outFile = $ARGV[5];

#Adding timestamp to jobID
my $timestamp = localtime(time);
$jobID = $jobID . $timestamp;

#Filling species hash, fill with all species taht COULD be asked to be found, should not change unless program overhaul
my %species = ();
$species{"Dana"} = "ananassae";
$species{"Dere"} = "erecta";
$species{"Dsec"} = "sechellia";
$species{"Dsim"} = "simulans";
$species{"Dyak"} = "yakuba";

#Fill hash with the species we WISH to search for on this run, if statements search for each type of species to be found, if found they are added to the hash
my %searchlist = ();
if($speciesString =~ /Ananassae/){ $searchlist{"Dana"} = 0; }
if($speciesString =~ /Erecta/){ $searchlist{"Dere"} = 0; }
if($speciesString =~ /Sechellia/){ $searchlist{"Dsec"} = 0; }
if($speciesString =~ /Simulans/){ $searchlist{"Dsim"} = 0; }
if($speciesString =~ /Yakuba/){ $searchlist{"Dyak"} = 0; }

#inFile should be formatted correctly, ie each name to be searched for on a new line, no header file
open(iFile, $inFile) || die "FAILED TO OPEN DATAFILE!\n";
my @geneNames = <iFile>;
close iFile;
chomp(@geneNames);

#Open file to save data to
open(OUTFILE, ">".$outFile) || die "FAILED TO OPEN DATAFILE!\n";

#Print header for file
print OUTFILE ">GENE_NAMES: @geneNames ; UP_FLANK_LENGTH: $upLen ; DOWN_FLANK_LENGTH: $downLen ; $timestamp\n##### START DATA FILE ##########################################################################################################\n"; 

#Query database with above parameters	
foreach $name (@geneNames){

	#Querying database
	my $sth = $dbh -> prepare("SELECT ID, to_species, to_id, to_name, target FROM gene WHERE name LIKE ?");
	$sth -> execute($name);
	
	print OUTFILE "##### GENE: $name DATA START ##################################################################################################\n"; 

 	#Fetch querys
	while(my @val = $sth->fetchrow_array()){
		
		# val[0] = ID
		# val[1] = to_species
		# val[2] = to_id
		# val[3] = to_name
		# val[4] = target

		#Setting curent species and target region
		my $curspecies = $val[1];
		my $curtarget = $val[4];
        	
		#Only process if the species exists in the searchlist hash which is genes we wish to find on this run
		if(exists $searchlist{$curspecies}){
			
			#Setting name of file to open and opening with BIOPERL
			my $path = "/home/max/Documents/galaxy/galaxy-dist/tools/myTools/";
			my $inputFile = $path . $species{$curspecies} . ".fa";
			print "inptuFile:$inputFile\n";
			my $seqio = Bio::SeqIO->new(-file => $inputFile, '-format' => 'Fasta');

			#Spitting target by spaces and setting values
			my @targetsplit = split(' ', $curtarget);	
			my $region = $targetsplit[0];
			my $num1 = $targetsplit[1];
			my $num2 = $targetsplit[2];
			my $direction = $targetsplit[3];
			
			#These variables hold the values from the search below
			my $string;
			my $id;
			my $desc;
			
			#cycle through the fasta file till we find the region we are looking for then use last to jump out of the loop
			while(my $seq = $seqio->next_seq){
				$id = $seq->id;
				if($id eq $region){
					$string = $seq->seq;
					$desc = $seq->desc;
					last;						
				}				
			}

			print OUTFILE "\n>FASTA DESC: $desc\n";
			my $to_name = substr $val[3], 5;
			
			#if upstream or downstream and print data to outfile
			if($direction eq "+"){
		
				my $upstream = substr $string, ($num1 - $upLen - 1), $upLen;
				my $downstream = substr $string, $num2, $downLen;
				print OUTFILE">FromName_$name"."_ToName_$to_name"."_ToSpecies_$curspecies" . "_" . ($num1 - $upLen - 1) . "_$upLen" . "_Upstream\n" . "$upstream\n";
				print OUTFILE">FromName_$name"."_ToName_$to_name" . "_ToSpecies_$curspecies" . "_" . $num2 . "_$downLen" . "_Downstream\n" . "$downstream\n";
				print OUTFILE"\n**********************************************************************************************************************\n";
				
			} else {

				my $downstream = substr $string, ($num1 - $downLen - 1), $downLen;
				my $upstream = substr $string, $num2, $upLen;
				print OUTFILE">FromName_$name"."_ToName_$to_name" . "_ToSpecies_$curspecies" . "_" . $num2 . "_$upLen" . "_Upstream\n" . "$upstream\n";
				print OUTFILE">FromName_$name"."_ToName_$to_name"."_ToSpecies_$curspecies" . "_" . ($num1 - $downLen - 1) . "_$downLen"."_Downstream\n"."$downstream\n";
				print OUTFILE"\n**********************************************************************************************************************\n";
				
			}
		
		}

	} #close fetch array from query loop
    print OUTFILE "\n##### END GENE  ###############################################################################################################################################################\n\n";
	$sth -> finish;
}

close OUTFILE;
