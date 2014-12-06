#! /usr/bin/perl

# This program parses all jobs.

use DBI; 
use DBD::mysql;

my $ds = "DBI:mysql:openMotif:localhost";
my $user = "root";
my $passwd = "egaga1234";
my $dbh = DBI->connect($ds, $user, $passwd) || die "Cannot open db";



print "Loading jobs into database. This may take a few moments.\n";

opendir(DIR, ".") or die "cannot open dir $dir: $!";
my @directory = readdir DIR;
closedir DIR;

my @newDirectories;

foreach my $row (@directory)
{
  chomp $row;
  if($row =~ "autoResults")
  {
    chomp $row;
    $directory = $row;
    my @jobFinder = split('_', $row);
    $job = "";
    for(my $i = 0; $i < scalar(@jobFinder); $i++)
    {
      if($i == scalar(@jobFinder) - 2)
      {
	$job .= $jobFinder[$i];
      }
      elsif($i < scalar(@jobFinder) - 2)
      {
        $job .= $jobFinder[$i] . "_";
      }
    }
  }
  my $flag = "bad";

  my $test = "SELECT COUNT(*) FROM Job WHERE job_id = ?";
  my $sth = $dbh->prepare($test);
  $sth->execute($job);

    my $ref = $sth->fetchrow_arrayref();

    if($ref->[0] != 0)
    {
      my $ans;
      print "The $job job is already in the database. Would you like to input this job under a different job ID?";
      $ans = <>;
      chomp $ans;
      if($ans eq "y" || $ans eq "yes" || $ans eq "Yes" || $ans eq "Y")
      {
        print "Please enter a new Job ID: ";
        my $newJob = <>; chomp $newJob;

        $sth -> execute($newJob);
        my $ref = $sth->fetchrow_arrayref();
        while($ref->[0] != 0)
        {
	  print "$ref->[0]\n";
	  print "That job is also already in the database. Try another job ID: ";
          $newJob = <>; chomp $newJob;
          $sth -> execute($newJob);
          my $ref = $sth->fetchrow_arrayref();
        }
        my $newDir = $newJob . "_autoResults";
        mkdir $newDir;
        opendir(DIR, $directory);
        my @fileSet = readdir DIR;
        closedir DIR;
        foreach $file (@fileSet)
        {
          my @finder = split($job, $file);
          
          if($file ne "." && $file ne "..")
          {
            chomp $file;
            open(INPUT, $directory . "/" . $file);
            open(OUTPUT, ">" . $newDir . "/" . $newJob . $finder[1]);
  
            my @file_contents = <INPUT>;
            foreach $line (@file_contents)
            {
              print OUTPUT "$line";
            }
          }
        }
        push(@directory, $newDir);
        push(@newDirectories, $newDir);
        next;
      }
      else
      {
        next;
      }
    }
    else
    {
  chomp $job;
  print "JOB: $job\n";

#if($row ne "." && $row ne ".." && $row ne "OpenMotif_DB_Loader.pl")
if($row =~ "autoResults")
{
  print "processing $row...\n";
chomp $row;

print "$directory / $job\n";
my $myfile = $directory . "/" . $job . "_wordClusters.pdf";

open MYFILE, $myfile or die "Cannot open file";
my $data;

while(<MYFILE>) {
  $data .= $_;
}

close MYFILE;

#Inserting into the Job table
print "Inserting into $job Job table...\n";
my $sql = "INSERT INTO Job VALUES(?,?)";
my $sth = $dbh->prepare($sql);
my $numrows = $sth->execute($job,$data);


#Now, I will parse the hm05_cluster_R_Cluster_Result.csv file
open(FILE, $directory . "/" . $job . "_R_Cluster_Result.csv");
@raw_clusters_file = <FILE>;
close(FILE);

my %cluster_hash;
$i = 0;
$len = scalar(@raw_clusters_file);



LOOP:
{
  while($i < $len)
  {
    if($raw_clusters_file[$i] =~ /Cluster/)
    {
      my $cluster_name = $raw_clusters_file[$i];
      chomp $cluster_name;
      $i++;
      while($raw_clusters_file[$i] ne "\n")
      {
        if($i >= $len)
        {
	  last LOOP;
        }
        my $word = $raw_clusters_file[$i];
        chomp $word;
        $cluster_hash{$word} = $cluster_name;
        $i++;
      }
    }
    $i++;
      if($raw_clusters_file[$i] eq "\n")
      {
        $i++;
      }
   }
}



#PARSING pfms and motif_score, storing into hashes.
open(FILE, $directory . "/" . $job . "_pfmLine.txt");
@raw_pfmLine = <FILE>;
close FILE;

my %motif_score;
my %pfm_a;
my %pfm_c;
my %pfm_g;
my %pfm_t;

for($i = 0; $i < scalar(@raw_pfmLine); $i++)
{
  if($i % 2 == 0)
  {
    my @parsed = split(',', $raw_pfmLine[$i]);
    my @word = split(':', $parsed[0]);
    my @motif_score = split(':', $parsed[1]);
  
    $word = $word[1]; chomp $word;
    my $motif_score = $motif_score[1]; chomp $motif_score;

    $motif_scores{$cluster_hash{$word}} = $motif_score;
    $last_word = $word;
  }
  else
  {
    my @parsed = split(':', $raw_pfmLine[$i]);
    $pfm_a{$cluster_hash{$word}} = $parsed[0];
    $pfm_c{$cluster_hash{$word}} = $parsed[1];
    $pfm_g{$cluster_hash{$word}} = $parsed[2];
    $pfm_t{$cluster_hash{$word}} = $parsed[3];
  }
}


#########PARSING_MOTIF_LOGOS##############

opendir(DIR, $directory) or die "cannot open dir $dir: $!";
my @file = readdir DIR;
closedir DIR;

my %motifLabelHash;
foreach $row (@file)
{
  if($row =~ ".png")
  {
  my @wordFinder = split("_", $row);
  my $word = $wordFinder[scalar(@wordFinder) - 2];
  $motifLabelHash{$cluster_hash{$word}} = $directory . "/" . $row;
  }
}

#CREATE A UNIQUE SET OF CLUSTERS.
my %unique_clusters;
foreach $key (keys %cluster_hash)
{
  $unique_clusters{$cluster_hash{$key}} = $cluster_hash{$key};
}

print "Inserting into $job Cluster table...\n";
foreach my $cluster (keys %unique_clusters)
{
  open MYFILE, $motifLabelHash{$cluster} or die "Cannot open file";
  my $data;

  while(<MYFILE>) {
    $data .= $_;
  }

  close MYFILE;

#Insertion into the cluster table
  my $sql = "INSERT INTO Cluster VALUES(?,?,?,?,?,?,?,?)";
  my $sth = $dbh->prepare($sql);
  my $numrows = $sth->execute($cluster,$pfm_a{$cluster},$pfm_c{$cluster},$pfm_g{$cluster},$pfm_t{$cluster},
			      $motif_scores{$cluster},$data,$job);
}

# Parsing out S, coverage from seq_cov.csv, storing in hashes

open(FILE, $directory . "/" . $job . "_seqCov.csv");
my @seqCov = <FILE>;
close FILE;


my %coverage_hash;
my %S_hash;

for($i = 2; $i < scalar(@seqCov); $i++)
{
  @parsed_seqCov = split(',', $seqCov[$i]);
  my $word = $parsed_seqCov[0]; chomp $word;
  $seq_word_count = $parsed_seqCov[1]; chomp $seq_word_count;
  $coverage = $parsed_seqCov[2]; chomp $coverage;
  $coverage_hash{$word} = $coverage;
  $S_hash{$word} = $seq_word_count;
}

# Parsing O (total occurrences) and E (expected count) from HammingClusters file, storing in hashes

open(FILE, $directory . "/" . $job . "_HammingClusters.csv");
@HammingArr = <FILE>;
close FILE;

my %occurrences_hash;
my %expected_count_hash;

foreach $row (@HammingArr)
{
  if($row =~ "seedword")
  {
    my @firstParse = split(',', $row);
    my @wordParse = split(':', $firstParse[0]);
    my $word = $wordParse[1];
    chomp $word;
    my @O_Parse = split(':', $firstParse[2]);
    chomp $O_Parse[1];
    $occurrences_hash{$word} = $O_Parse[1];
    my @E_Parse = split(': ', $firstParse[3]);
    chomp $E_Parse[1];
    $expected_count_hash{$word} = $E_Parse[1];
  }
}

#Parsing out avgE through Orders from hm05_clusterWords.csv, storing in hashes
open(FILE, $directory . "/" . $job . "_clusterWords.csv");
my @clusterWords = <FILE>;
close FILE;

my %S;
my %avgE;
my %avgSln;
my %avgZ;
my %avgPval;
my %stdE;
my %stdSln;
my %stdZ;
my %stdPval;
my %orders;

$i = 0;
while($i < scalar(@clusterWords))
{
  if($clusterWords[$i] =~ "avgE")
  {
    $i++;
    while($clusterWords[$i] ne "\n" && $i < scalar(@clusterWords))
    {
      my @parsed_line = split(',', $clusterWords[$i]);
      $word = $parsed_line[0];
      $S{$word} = $parsed_line[1];
      $avgE{$word} = $parsed_line[3];
      $avgSln{$word} = $parsed_line[4];
      $avgZ{$word} = $parsed_line[5];
      $avgPval{$word} = $parsed_line[6];
      $stdE{$word} = $parsed_line[7];
      $stdSln{$word} = $parsed_line[8];
      $stdZ{$word} = $parsed_line[9];
      $stdPval{$word} = $parsed_line[10];
      $orders{$word} = $parsed_line[11];
      $i++;
    }
  }
  $i++;
}


#Parsing out word counts, storing in a hash
open(FILE, $directory . "/" . $job . "_wordPos.csv");
my @wordPositions = <FILE>;
close FILE;

my %wordPosHash;
my @sequences;
my @seedSequences = split(',', @wordPositions[0]);

for(my $i = 1; $i < scalar(@seedSequences) - 1; $i++)
{
  my @temp = split('>', $seedSequences[$i]);
  my @temp2 = split(':', $temp[1]);
  chomp $temp2[0];
  push(@sequences, $temp2[0]);
}

for(my $index = 0; $index < scalar(@wordPositions); $index++)
{
  if($wordPositions[$index] =~ /seedword/)
  {
    for(my $i = 0; $i < scalar(@sequences); $i++)
    {
      $seq = $sequences[$i];
      chomp $seq;
      my @line = split(',', $wordPositions[$index]);
      my @wordFinderArr = split(':', $line[0]);
      my $word = $wordFinderArr[1]; chomp $word;
      my $theKey = $word . "," . $seq;
      @seqPositions = split(':', $line[$i+1]);

      if($seqPositions[0] ne "none")
      {
        my @shiftedPos;
	for(my $h = 1; $h < scalar(@seqPositions); $h++)
	{
	  push(@shiftedPos, $seqPositions[$h]);
	}
	$wordPosHash{$theKey} = \@shiftedPos;
      }
    }
  }
}


#Parsing out motif positions, storing into a hash
my %motifPosHash;
my %motifCountHash;

open(FILE, $directory . "/" . $job . "_motifPos.csv");
my @motifPositions = <FILE>;
close FILE;

for(my $i = 1; $i < scalar(@motifPositions); $i++)
{
  my @splitLine = split(',', $motifPositions[$i]);
  for(my $i = 0; $i < scalar(@sequences); $i++)
  {
    if($splitLine[$i+4] ne "none")
    {
      my @positions = split(':', $splitLine[$i+4]);
      my $motifCount = scalar(@positions);
      $motifPosHash{$cluster_hash{$splitLine[0]} . "," . $sequences[$i]} = \@positions;
      $motifCountHash{$cluster_hash{$splitLine[0]} . "," . $sequences[$i]} = $motifCount;
    }
    else
    {
      $motifCountHash{$cluster_hash{$splitLine[0]} . "," . $sequences[$i]} = 0;
      my $key = $cluster_hash{$splitLine[0]} . "," . $sequences[$i];
    }
  }
}

open(FILE, $directory . "/" . $job . "_wordCount.csv");
my @wordCounts = <FILE>;
close FILE;

my %wordCountHash;
for(my $index = 1; $index < scalar(@wordCounts); $index++)
{
  if($wordCounts[$index] =~ "seedword")
  {
    my @splitCounts = split(",", $wordCounts[$index]);
    my @wordFindArr = split(":", $splitCounts[0]);
    my $word = $wordFindArr[1];
    for(my $i = 0; $i < scalar(@sequences); $i++)
    {
      $wordCountHash{$word . "," . $sequences[$i]} = $splitCounts[$i+1];
    }
  }
}

#Inserting into Sequence table
print "Inserting into $job Sequence table...\n";
foreach $row (@sequences)
{
  my $sql = "INSERT INTO Sequence VALUES(?,?)";
  my $sth = $dbh->prepare($sql);
  my $numrows = $sth->execute($row, $job);
}

#Inserting into Word table
print "Inserting into $job Word table...\n";
foreach $key (keys %cluster_hash)
{
  my $job_id = $job;
  my $sql = "INSERT INTO Word VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
  my $sth = $dbh->prepare($sql);
  my $insert = $sth->execute($key,$cluster_hash{$key},$S{$key},
	   		     $coverage_hash{$key},$occurrences_hash{$key},
			     $avgE{$key},$avgSln{$key},$avgZ{$key},$avgPval{$key},
			     $stdE{$key},$stdSln{$key},$stdZ{$key},$stdPval{$key},
			     $orders{$key},$expected_count_hash{$key},$job_id);
}  

#Inserting into Word Position table
print "Inserting into $job Word_pos table...\n";
foreach my $key (keys %wordPosHash)
{
  foreach my $position (@{$wordPosHash{$key}})
  {
    my $job_id = $job;
    my @finder = split(",", $key);
    my $word = $finder[0];
    my $seq = $finder[1];
    
    my $sql = "INSERT INTO Word_pos VALUES(?,?,?,?)";
    my $sth = $dbh->prepare($sql);
    my $insert = $sth->execute($word,$seq,$position,$job_id);
  }
}

#Inserting into Motif Position table
print "Inserting into $job Motif_pos table...\n";
foreach my $key (keys %motifPosHash)
{
  foreach my $position (@{$motifPosHash{$key}})
  {
    my $job_id = $job;
    my @finder = split(",", $key);
    my $cluster_name = $finder[0];
    my $seq = $finder[1];
    
    my $sql = "INSERT INTO Motif_pos VALUES(?,?,?,?)";
    my $sth = $dbh->prepare($sql);
    my $insert = $sth->execute($cluster_name,$seq,$position,$job_id);
  }
}

print "Inserting into $job Motif_count table...\n";
foreach $key (keys %motifCountHash)
{
    my $job_id = $job;
    my @finder = split(",", $key);
    my $cluster_name = $finder[0];
    my $seq = $finder[1];

    my $sql = "INSERT INTO Motif_count VALUES(?,?,?,?)";
    my $sth = $dbh->prepare($sql);
    my $insert = $sth->execute($cluster_name,$seq,$motifCountHash{$key},$job);
}
  
#Inserting into Word Count table
print "Inserting into $job Word_count table...\n";
foreach my $key (keys %wordCountHash)
{
  {
     my @finder = split(",", $key);
     my $word = $finder[0];
     my $seq = $finder[1];
     my $count = $wordCountHash{$key};
     my $job_id = $job;

    my $sql = "INSERT INTO Word_count VALUES(?,?,?,?)";
    my $sth = $dbh->prepare($sql);
    my $insert = $sth->execute($word,$seq,$count,$job_id);
  }
}
}
}
}

print "Files uploaded.\n";
