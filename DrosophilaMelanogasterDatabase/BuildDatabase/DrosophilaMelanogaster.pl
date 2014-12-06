#! /usr/bin/perl
######################################################################################################
#
#	Program: 	DrosophilaMelanogaster.pl
#	Version:	1.0
#
#	Author:		Max Rego
#	Email:		mr255509@ohio.edu
#	Project:	Drosophila Melanogaster Orthologous Genes Database
#	
#	Description:	This program parses a datafile from Flybase and extracts useful information
#			for a database to find orthogous genes
#
#	Date:		3-14-2013
#
######################################################################################################
#
################  USE SECTION   ######################################################################

use DBI;
use DBD::mysql;

################   ORTH CLASS START  #################################################################

package orth;
sub new
{
    my $class = shift;
    my $self = {
        _NAME  		=> shift,
        _to_species 	=> shift,
        _to_id      	=> shift,
        _to_name    	=> shift,
        _target		=> shift,
    };    
   
    bless $self, $class;
    return $self;
}

sub getNAME {
    my( $self ) = @_;
    return $self->{_NAME};
}

sub getTo_species {
    my( $self ) = @_;
    return $self->{_to_species};
}

sub getTo_id {
    my( $self ) = @_;
    return $self->{_to_id};
}

sub getTo_name {
    my( $self ) = @_;
    return $self->{_to_name};
}

sub getTarget {
    my( $self ) = @_;
    return $self->{_target};
}

############################################################################################################

#########    MAIN    #######################################################################################

#Connecting to Database
my $ds = "DBI:mysql:dmel:localhost";
my $user = "";
my $passwd = "";
my $dbh = DBI -> connect($ds, $user, $passwd) || die "CANT CONNECT!\n";

#Opening File - file array is @file, $filesize is the number of lines in the file
open(GENES, "dmel.gff") || die "FAILED TO OPEN DATAFILE!\n";
my @file = <GENES>;
close GENES; 
my $filesize = @file;
print "File opened!\n";

#Global Varibales
my %sequence_regions = ();	#holds the sequence regions from the beginning of the .gff file
my $i = 0;			#holds the number of the current working line
my %genes = (); 		#holds the gene events
my $orthcount = 0;		#holds the number of orthologous_to's triggered
my %species = ();

#Filling species hash
$species{"Dana"} = "ananassae";
$species{"Dere"} = "erecta";
$species{"Dsec"} = "sechellia";
$species{"Dsim"} = "simulans";
$species{"Dyak"} = "yakuba";

#Parse File By Each Line
foreach $line (@file){

	if ( $line =~ /##/ ) {	#checking for sequence regions and comment lines

		if($line =~ /##sequence-region/){	#checking if the section is a sequence-region
			
			my @sequence_info = split(' ', $line);
			$sequence_regions{$sequence_info[1]} = $sequence_info[3];
		
		} elsif ($line =~ /##FASTA/){	

			#Jump out of loop if ##FASTA appears because thats when the fasta files start
			last;

		} else {

			#other comment region

		}
	
	} elsif ( ($line =~ /orthologous_to/) ){  #contains orthogous_to which is the section of data we are looking for

		my @linesplit = split('\t', $line);	#splitting the current line by tabs, always should be 9 pieces
		my $linenum = @linesplit;		#number of pieces in the line
		$orthcount++;				#adding 1 to orthogous_to count

		#temp variables
		my $ID;
		my $NAME;
		my $to_species;
		my $to_id;
		my $to_name;
		my $target = "NONE";	#Assigned to NONE in case the orthologous section doesnt have a Target keyword in its info section

		#Beginning search for info section
		#I didnt want to hard code which section in @linesplit holds the info section in case the data file is ever changed.
		#Currently linesplit[8] should always hold the info section but i have my algorithm searching for it just in case.
		#After finding which $split of @linesplit holds the info section I split the info section by ; which is how its divided up
		#Each $sec is run through an if elsif switch below to find if it conatins a keyword that we are looking for
		#If found the part after the keyword is assigned to a temp variable and chomped, then assigned to the variable to hold that value

		foreach $split (@linesplit){
			if($split =~ /ID=/){
				my @infosplit = split(';', $split);
				foreach $sec (@infosplit){
					if($sec =~ /ID=/){

						my $temp = $';
						chomp($temp);
						$ID = $temp;

					} elsif ($sec =~ /Name=/){

						my $temp = $';
						chomp($temp);
						$NAME = $temp;

					} elsif( $sec =~ /to_species=/){

						my $temp = $';
						chomp($temp);
						$to_species = $temp;

					} elsif ($sec =~ /to_id=/){

						my $temp = $';
						chomp($temp);
						$to_id = $temp;

					} elsif ($sec =~ /to_name=/){

						my $temp = $';
						chomp($temp);
						$to_name = $temp;

					} elsif ($sec =~ /Target=/){

						my $temp = $';
						chomp ($temp);
						$target = $temp;

					}

				}
			}
		} 

		#Storing Current Orthologus_to Line in a Class, Then Storing in Genes Hash with its ID as its Unquie Identifier
		$object = new orth($NAME, $to_species, $to_id, $to_name, $target);
		$genes{$ID} = $object;

	}

	$i++;

}

#Loop Through Genes Hash to Add to the Database
print "Starting DB addition\nOrthcount = $orthcount\n";
while (($key, $value) = each(%genes)){

	#assigning temp variables to hold each data value for each key
	my $ID = $key;
	my $NAME = $value->getNAME();
	my $to_species = $value->getTo_species();
	my $to_id = $value->getTo_id();
	my $to_name = $value->getTo_name();
	my $target = $value->getTarget();
    
	#inserting into database
	my $sth = $dbh -> prepare("INSERT INTO gene(ID,name,to_species,to_id,to_name,target) values(?,?,?,?,?,?)");
	$sth -> execute ($ID, $NAME, $to_species, $to_id, $to_name, $target);
	$sth -> finish;

}

print "Database Loaded\n";

#undef global variables and disconnecting from database
undef %genes;
undef %sequence_regions;
undef $i;
$dbh->disconnect;




