#################################################################################################
#
#	Program:	query_menu.pl
#
#	Author:		Max Rego
#	Email:		mr255509@ohio.edu
#	
#	Description:	This program is the menu for querying the database created to store the
#			output from OpenMotif and generate a PDF report containing the results
#
#	Note:		Must have PDF::API2 installed from CPAN
#
#	Date:		12-3-12
#
#################################################################################################

#! /usr/bin/perl

use PDF::API2;
use DBI;
use DBD::mysql;

#Insert database login info here 
my $ds = "DBI:mysql:    :localhost";
my $user = "   ";
my $passwd = "   ";
my $dbh = DBI->connect($ds, $user, $passwd) || die "Cannot open Database";

my $outer_flag = 1;		#controls the query menu loop; primed to 1; if anything else then loop stops 
my $switch_cont = 0; 		#controls the if statement switch inside outer loop 
my $query_string;		#holds the query and parameters entered by the user; format {query_number;} | {query_number, parameter_x;} | {query_number ,parameter_x, parameter_y;} 
				#note that each valid query will be added to the string, parameters may not be needed, example in query 2

#start query menu; querys will be added to a string to be executed, user enters -1 to finish 
print "\n* * * * * * * * * * * * * * * *  * * * * * * * * QUERY MENU * * * * * * * * * * * * * * * * * * * * * * * * * * * \n";

while($outer_flag == 1){

	# Prints querys already in the query_string, will only print querys with 0 1 or 2 parameters, if more is needed add another condition in the loop below. 
	print "\n* * * So far these querys have been added to the report:\n";
	my $tprint_string = $query_string;
	my @tprint_arr = split(';', $tprint_string);
	foreach $tquery (@tprint_arr){
		my @tquery_bits = split(',', $tquery);
		my $tpos = 1;		
		foreach $bit (@tquery_bits){
			if($tpos == 1){
				print "Query number $bit ";
				$tpos++;
			} elsif ($tpos = 2){
				print "with Parameter $bit ";
				$tpos++;
			} else {
				print "and with parameter $bit";
				$tpos++;
			}
		}
		print "\n";
	}	

	#Displays the menu to user
	print "\n* * * Please choose an option to add to the report:\n";
	print " 1. Enter a specific word (string) and check if it exists in the database.\n";
	print " 2. Return the word with the largest number of occurances.\n";
	print " 3. Return all words that occur between a pair of coordinates.\n";
	print " 4. Return the input sequences from which a given word occurs in.\n";
	print " 5. Return the motif with the largest number of occurances.\n";
	print " 6. Return the 5 highest scoring motifs.\n";
	print " 7. Return all words that belong to a certain motif.\n";
	print " 8. Return the motifs with the highest sequence coverage.\n";
	print "-1. Generate Report and Exit.\n";
	print ">> ";
	$switch_cont = <>;
	print "\n";			

	#Begin if switch to add to query_string	
	#To add a new query simply add another elsif statement at the end before final else -> elsif($switch_cont == X)
	#Most code can be reused for new querys with 0 1 or 2 parameters; 
	#*******Safety feature not included, if a user enters a ; in the statement it WILL cause trouble when querying, should implement a check for user entered ;  
	if($switch_cont == 1){
		print "You have entered query 1\n";
		print "Please enter the word you wish to check or -1 to cancel: ";		
		my $temp_word = <>;
		if($temp_word != -1){
			chomp($temp_word);
			$query_string = $query_string . "1," . $temp_word . ";";
			print "Successfully added to report!....Relaoding Menu\n"; 
		} else {
			print "Query Cancelled! Nothing was added to the report...Reloading Menu\n";
		}
	} elsif ($switch_cont == 2){
		if($query_string =~ /;2;/){	#checking if query 2 has been added to the report already
			print "Sorry Query 2 has been added to the report already...Reloading Menu\n";
		} else {
			print "You have entered query 2\n";
			print "Successfully added to the report!...Reloading Menu\n";
			$query_string = $query_string . "2;";	
		}
	} elsif ($switch_cont == 3){
		print "You have entered query 3\n";
		my $temp_coord_flag = 1;	#primed to 1, only changes if user enters valid coordinates or -1 to cancel addition to report
		my $temp_lcoord;		#holds the lower end coordinate
		my $temp_ucoord;		#holds the upper end coordinate	
		while($temp_coord_flag == 1){		
			print "Please enter the lower end coordinate or -1 to cancel: ";
			$temp_lcoord = <>;
			print "\nPlease enter an upper end coordinate or -1 to cancel: ";
			$temp_ucoord = <>;
			if ($temp_lcoord == -1 || $temp_ucoord == -1) {		#checking if the user cancelled the addition
				print "\nQuery Cancelled! Nothing was added to the report...Reloading Menu\n";
				$temp_coord_flag = 0;
			} elsif	($temp_lcoord >= $temp_ucoord){			#checking if coordinates are valid
				print "\nInvalid coordinates (make sure lower end coordinate is not the same as or larger than the upper end coordinate)...try again\n";
			} else {						#if coords are valid and no cancellation was requested exit loop
				$temp_coord_flag = 0;
			}
		}
		if($temp_lcoord != -1 && $temp_ucoord != -1){	#check for cancellation, if no cancellation then add query to query_string
			chomp($temp_lcoord);
			chomp($temp_ucoord);			
			$query_string = $query_string . "3," . $temp_lcoord . "," . $temp_ucoord . ";";
			print "Successfully added to the report!...Reloading Menu\n";		
		} else {
			print "\nQuery Cancelled! Nothing was added to the report...Reloading Menu\n";
		}			
	} elsif ($switch_cont == 4){
		print "You have entered query 4\n";
		print "Please enter the Word to find which sequences it occurs in or -1 to cancel: ";
		my $temp_word = <>;
		chomp($temp_word);
		if($temp_word != -1){	#Checks for cancellation, if no cancellation then add query to query_string
			$query_string = $query_string . "4," . $temp_word . ";";			
			print "Successfully added to the report!...Reloading Menu\n";
		} else {
			print "\nQuery Cancelled! Nothing was added to the report...Reloading Menu\n";
		}
	} elsif ($switch_cont == 5){
		if($query_string =~ /;5;/){     #checking if query 5 has been added to the report already
                        print "Sorry Query 5 has been added to the report already...Reloading Menu\n";
                } else {
                        print "You have entered query 5\n";
                        print "Successfully added to the report!...Reloading Menu\n";
                        $query_string = $query_string . "5;"; 
                }
	
	} elsif ($switch_cont == 6){
                if($query_string =~ /;6;/){     #checking if query 6 has been added to the report already
                        print "Sorry Query 6 has been added to the report already...Reloading Menu\n";
                } else {        
                        print "You have entered query 6\n";                     
                        print "Successfully added to the report!...Reloading Menu\n";
                        $query_string = $query_string . "6;";
                } 

	} elsif ($switch_cont == 7){
		print "You have entered query 7\n";
                print "Please enter the motif you wish to find words for or -1 to cancel: ";
                my $temp_word = <>;
                if($temp_word != -1){
                        chomp($temp_word); 
                        $query_string = $query_string . "7," . $temp_word . ";";
                        print "Successfully added to report!....Relaoding Menu\n";
                } else {        
                        print "Query Cancelled! Nothing was added to the report...Reloading Menu\n";
                }

	} elsif ($switch_cont == 8){
		if($query_string =~ /;8;/){     #checking if query 8 has been added to the report already
                        print "Sorry Query 8 has been added to the report already...Reloading Menu\n"; 
                } else {
                        print "You have entered query 8\n";
                        print "Successfully added to the report!...Reloading Menu\n";
                        $query_string = $query_string . "8;";
                }

	} elsif ($switch_cont == -1){
		print "Generating report now...\n"; #exit loop and generate report
		$outer_flag = 0;

	} else {
		print "Invalid entry!...No query was executed!...Reloading menu\n";
	}
}

#declaring a new PDF
$pdf = PDF::API2->new();

	# NOTES ON PDF FILES
	# -output to PDF does not text wrap, if text string to print to pdf is too long it will run off the page
	#  Not implemented but to fix keep a variable to watch the length of the string and make sure its not too long
	# -module can be buggy, its not documented extremely well, to keep anything from going wrong I copied 
	#  things at each print statement that may not be needed every time, see below example 
	#
	# -Here is an example format for printing a text line
	# 
	#  $page = $pdf->page();		//adds a new page to the pdf
        #  $page->mediabox('Letter');		//sets page format limits
        #  $font = $pdf->corefont('Helvetica-Bold');	//set font, can add new fonts with different function, look at PDF::API2 documentation
        #  $text = $page->text();		//prepare a text line
        #  $text->font($font, 40);		//set font and size:  $text->font(font, size)	
        #  $text->translate(70, $cursor);	//set location to print at:   $text->translate(X location, Y location) SEE NOTE BELOW ON FORMAT OF LOCATIONS 
        #  $text->text($report_name);		//prints the string to the pdf at above location, $text->test($string) | $text->test('Hello world') 
        #  $cursor = $cursor - 20;		//set new cursor
	#
	#  FORMAT OF LOCATIONS
	#  -X location -> distance from left side of page, I found location 70 is a good place to start
	#  -Y location -> distance from bottom of page, I found 700 is a good place to start, -20 for each text line, -220 for each image	
	#
	#  PICTURE FORMAT
	#  -loading pictures from the database is done by downloading the string and creating a .png file in the current directory, then adding the image
	#   to the PDF, then deleting the .png 
	#
	#  -Heres an example of printing an image		
	#
	#  my $temp_file = "temp.png";		//variable for temp file	
	#  open(MYOUTFILE, ">$temp_file");	//create a new file with the file variable   
        #  print MYOUTFILE $temp_image;		//print the string fetched from the database
        #  close(MYOUTFILE);			//close it
        #  my $gfx = $page -> gfx;			//starts a graphics image on the page
        #  my $image = $pdf -> image_png($temp_file);	//fetching the string from the temp file made above
        #  $gfx -> image($image, 100, 460, 400, 200);	//format and print image to PDF; $gfx->image($image, X_Coord, Y_Coord, X_Image_Size, Y_Image_Size)
	#						//X_Coord = location from the left side of page; Y_Coord = location from the bottom of the page
	#						//X_Image_Size = width of the image across; Y_Image_Size = height of the image
	#						//i found that (400, 200) is a good size for the images from OpenMotif
	# unlink($temp_file);				//delete the .png created above					
	

	#Naming the report
	print "Pleass give a name for the Report (do not include an extention): ";
	my $report_name = <>;
	chomp($report_name);
	
	#printing header to PDF
	my $cursor = 640;
        $page = $pdf->page();
        $page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
	$text = $page->text();
        $text->font($font, 40);
       	$text->translate(70, $cursor);
	$text->text($report_name);
	$cursor = $cursor - 20;
	
	#Add OU Picture
	my $gfx = $page -> gfx;
        my $image = $pdf -> image_png('testw.png');
        $gfx -> image($image, 80, 150, 450, 450);

	#Add timestamp
	@months = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
 	@weekDays = ("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun");
 	($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime(); 
	$year = 1900 + $yearOffset;
 	$theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";	
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 20);
        $text->translate(70, 75);
        $text->text($theTime);
        $cursor = $cursor - 20;

	#print the query table
	$cursor = 700;
	$page = $pdf->page();
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 30);
        $text->translate(70, $cursor);
        $text->text('Query Table');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('1. Enter a specific word (string) and check if it exists in the database.');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('2. Return the word with the largest number of occurances.');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('3. Return all words that occur between a pair of coordinates.');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('4. Return the input sequences from which a given word occurs in.');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('5. Return the motif with the largest number of occurances.');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('6. Return the 5 highest scoring motifs.');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('7. Return all words that belong to a certain motif.');
	$cursor = $cursor - 20;
	$page->mediabox('Letter');
        $font = $pdf->corefont('Helvetica-Bold');
        $text = $page->text();
        $text->font($font, 12);
        $text->translate(70, $cursor);
        $text->text('8. Return the motifs with the highest sequence coverage.');
	$cursor = $cursor - 20;

# Starting Querys
#Adding new querys shouldnt take too long, a lot of the code can be copied and pasted for quick addition of a query
my @fquery_arr = split(';', $query_string);
	foreach $fquery (@fquery_arr){
		my @fquery_bits = split(',', $fquery);
		if($fquery_bits[0] == 1){
			
			#querying database
			my $tpara = $fquery_bits[1];			
			my $sth = $dbh -> prepare("SELECT word, position, seq, job_id FROM Word_pos WHERE word LIKE ? ORDER BY job_id ASC"); 
        		$sth -> execute($tpara);
			
			#printing header to PDF
			my $cursor = 700;
                        $page = $pdf->page();
                        $page->mediabox('Letter');
                        $font = $pdf->corefont('Helvetica-Bold');
			$text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
			$text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 1 * * * * * * * * * * * * * * * * * * * * * * * * *');
			$cursor = $cursor - 20;
			$text = $page->text(); 
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
			my $print_para = "* * * * Using Parameter $tpara";
			$text->text($print_para); 
                        $cursor = $cursor - 20;
			
			#printing query data to PDF
			while(my @val = $sth->fetchrow_array()){ 
				my $temp_result = "Word: $val[0] | Position: $val[1] | Seq: $val[2] | Job ID: $val[3] \n";
				if($cursor > 49){
					$text = $page->text();
                                	$text->font($font, 12);
                                	$text->translate(70, $cursor);
                                	$text->text($temp_result);
					$cursor = $cursor - 20;
				} else {
					$cursor = 700;
					$page = $pdf->page();
                        		$page->mediabox('Letter');
                        		$font = $pdf->corefont('Helvetica-Bold');
					$text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result);
                                        $cursor = $cursor - 20;
				}
			}
		
			#printing closing line of query
			$text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');		
        		$sth -> finish;
	
		} elsif ($fquery_bits[0] == 2){
			#querying database			
			my $sth = $dbh -> prepare("SELECT title, occurrences, job_id FROM Word ORDER BY occurrences DESC"); 
        		$sth -> execute();
			
			#printing header of query to PDF
			my $cursor = 700;
                        $page = $pdf->page();
                        $page->mediabox('Letter');
                        $font = $pdf->corefont('Helvetica-Bold');
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 2 * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $cursor = $cursor - 20;
                        $text = $page->text(); 
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        my $print_para = "* * * * Word with most number of occurrences";
                        $text->text($print_para);
                        $cursor = $cursor - 20;

			#printing query data to PDF
			my @val = $sth->fetchrow_array();
				my $temp_result = "Word: $val[0] | Occurrences: $val[1] | job_id: $val[2] \n";
				$text = $page->text(); 
                        	$text->font($font, 12);
                        	$text->translate(70, $cursor); 
                        	$text->text($temp_result);
                        	$cursor = $cursor - 20;
        	
			#printing closing line of query
			$text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $sth -> finish;

		} elsif ($fquery_bits[0] == 3){
			
			#querying database
			my $tpara_x = $fquery_bits[1];
			my $tpara_y = $fquery_bits[2];		
			my $sth = $dbh -> prepare("SELECT word, position, seq, job_id FROM Word_pos WHERE position > ? AND position < ? ORDER BY position ASC"); 
        		$sth -> execute($tpara_x, $tpara_y);
				
			#printing header to PDF
	                my $cursor = 700;
       	  	        $page = $pdf->page(); 
       	     	        $page->mediabox('Letter');
                       	$font = $pdf->corefont('Helvetica-Bold');
                       	$text = $page->text();
                       	$text->font($font, 12);
                       	$text->translate(70, $cursor);
                       	$text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 3 * * * * * * * * * * * * * * * * * * * * * * * * *');
                       	$cursor = $cursor - 20;
                       	$text = $page->text(); 
                       	$text->font($font, 12);
                       	$text->translate(70, $cursor); 
                        my $print_para = "* * * * Using Parameters  X = $tpara_x and Y = $tpara_y";
                       	$text->text($print_para);
                 	$cursor = $cursor - 20;

			#printing query data to PDF
        		while(my @val = $sth->fetchrow_array()){ 
				my $temp_result = "Word: $val[0] | Position: $val[1] | Seq: $val[2] | Job ID: $val[3] \n";
				if($cursor > 49){
                                        $text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result); 
                                        $cursor = $cursor - 20;
                                } else { 
                                        $cursor = 700;
                                        $page = $pdf->page();
                                        $page->mediabox('Letter');
                                        $font = $pdf->corefont('Helvetica-Bold');
                                        $text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result);
                                        $cursor = $cursor - 20;
                                }
                        }
			
			#printing closing line of query
                        $text = $page->text(); 
                        $text->font($font, 12);
                        $text->translate(70, $cursor); 
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $sth -> finish;
		
		} elsif ($fquery_bits[0] == 4){

			#querying database
			my $tpara = $fquery_bits[1];			
			my $sth = $dbh -> prepare("SELECT word, seq, job_id FROM Word_count WHERE word LIKE ? ORDER BY seq ASC"); 
        		$sth -> execute($tpara);
			
			#printing header to PDF
                        my $cursor = 700;
                        $page = $pdf->page(); 
                        $page->mediabox('Letter');
                        $font = $pdf->corefont('Helvetica-Bold');
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 4 * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $cursor = $cursor - 20;
                        $text = $page->text(); 
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        my $print_para = "* * * * Using Parameters $tpara"; 
                        $text->text($print_para);
                        $cursor = $cursor - 20;

			#printing query data to PDF
        		while(my @val = $sth->fetchrow_array()){ 
				my $temp_result = "Word: $val[0] | Seq: $val[1] | Job ID: $val[2] \n";
				if($cursor > 49){
                                        $text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result);
                                        $cursor = $cursor - 20;
                                } else {
                                        $cursor = 700;
                                        $page = $pdf->page();
                                        $page->mediabox('Letter');
                                        $font = $pdf->corefont('Helvetica-Bold');
                                        $text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result);
                                        $cursor = $cursor - 20;
                                }
                        }

			#printing closing line to PDF
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $sth -> finish;

		 } elsif ($fquery_bits[0] == 5){

			#querying the database for the cluster
			my $cluster;	#holds best motif cluster
			my $count; 	#holds the number of occurances
			my $sth = $dbh -> prepare("SELECT cluster_name, count FROM Motif_count ORDER BY count DESC LIMIT 1");
                        $sth -> execute();
				#fetching cluster
                       		my @val = $sth->fetchrow_array();
				$cluster = $val[0];
				$count = $val[1];
                       	$sth -> finish;
			
			#querying the database for the motif label
			my $sth = $dbh -> prepare("SELECT motif_label FROM Cluster WHERE cluster_name LIKE ?");
			$sth -> execute($cluster);
				my @val = $sth->fetchrow_array();
				#saving image to string
				my $temp_image = $val[0];	             
				my $temp_file = "temp.png";	
			
			#saving image to a temporary .png, will be deleted later
			open(MYOUTFILE, ">$temp_file");
			print MYOUTFILE $temp_image; 
			close(MYOUTFILE); 
	
			#printing header to PDF
			$page = $pdf->page();
			$page->mediabox('Letter');
			$font = $pdf->corefont('Helvetica-Bold');
			$text = $page->text();
			$text->font($font, 12);
			$text->translate(70, 700);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 5 * * * * * * * * * * * * * * * * * * * * * * * * *');
			$text = $page->text(); 
                        $text->font($font, 12);
                        $text->translate(70, 680); 
			my $print_para = "* * * * The highest number of occurances is $count";
                        $text->text($print_para);

			#printing image to PDF
			my $gfx = $page -> gfx;
			my $image = $pdf -> image_png($temp_file);
			$gfx -> image($image, 100, 460, 400, 200);
				
			#printing closing line to PDF
			$text = $page->text();
                        $text->font($font, 12); 
                        $text->translate(70, 440);
			$text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'); 
				
			#deleting temporary .png
			unlink($temp_file);
                        $sth -> finish;
		} elsif ($fquery_bits[0] == 6){
			
			#querying database
                        my $sth = $dbh -> prepare("SELECT cluster_name, motif_score ,motif_label FROM Cluster ORDER BY motif_score DESC LIMIT 5");
                        $sth -> execute();
                                
                        #printing header to PDF 
                        my $cursor = 700;
                        $page = $pdf->page();
                        $page->mediabox('Letter');
                        $font = $pdf->corefont('Helvetica-Bold');
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 6 * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $cursor = $cursor - 20;
                        
			#printing query data to PDF
                        while(my @val = $sth->fetchrow_array()){
                                my $temp_result = "Cluster Name: $val[0] | Motif Score: $val[1]\n";
                                
				if($cursor > 49){ 
                                        $text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result); 
                                        $cursor = $cursor - 20;
                                } else {
                                        $cursor = 700; 
                                        $page = $pdf->page();
                                        $page->mediabox('Letter');
                                        $font = $pdf->corefont('Helvetica-Bold');
                                        $text = $page->text(); 
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result);
                                        $cursor = $cursor - 20;
                                }
	
				if($cursor > 259){
					#setting cursor to correct location
					$cursor = $cursor - 200;
					my $temp_image = $val[2];
                                	my $temp_file = "temp.png";
					
					#saving image to a temporary .png, will be deleted later
                        		open(MYOUTFILE, ">$temp_file");
                        		print MYOUTFILE $temp_image;
                        		close(MYOUTFILE);
					
					#printing image to PDF
                       			my $gfx = $page -> gfx;
                        		my $image = $pdf -> image_png($temp_file);
                        		$gfx -> image($image, 100, $cursor, 400, 200);
					$cursor = $cursor - 20;
		
					#deleting temporary .png
		                        unlink($temp_file);
				} else {
                        		$cursor = 700;
                                        $page = $pdf->page();

					#setting cursor to correct location
                                        $cursor = $cursor - 200;
                                        my $temp_image = $val[2];
                                        my $temp_file = "temp.png";
                        
                                        #saving image to a temporary .png, will be deleted later
                                        open(MYOUTFILE, ">$temp_file");
                                        print MYOUTFILE $temp_image;
                                        close(MYOUTFILE);
                                        
                                        #printing image to PDF 
                                        my $gfx = $page -> gfx;
                                        my $image = $pdf -> image_png($temp_file);
                                        $gfx -> image($image, 100, $cursor, 400, 200);
                                        $cursor = $cursor - 20;
                                        
                                        #deleting temporary .png
                                        unlink($temp_file);
				}
			}
                        
                        #printing closing line of query
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $sth -> finish;
				
		} elsif ($fquery_bits[0] == 7){
			
			#querying database
                        my $tpara = $fquery_bits[1];
                        my $sth = $dbh -> prepare("SELECT title, cluster_name, occurrences, job_id FROM Word WHERE cluster_name LIKE ? ORDER BY job_id ASC");
                        $sth -> execute($tpara);
                 
                        #printing header to PDF 
                        my $cursor = 700;
                        $page = $pdf->page();
                        $page->mediabox('Letter');
                        $font = $pdf->corefont('Helvetica-Bold');
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 7 * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $cursor = $cursor - 20;
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor); 
                        my $print_para = "* * * * Using Parameter $tpara";
                        $text->text($print_para);
                        $cursor = $cursor - 20;
			
			#printing query data to PDF
                        while(my @val = $sth->fetchrow_array()){
                                my $temp_result = "Title: $val[0] | Cluster Name: $val[1] | Occurrences: $val[2] | Job ID: $val[3] \n";
                                if($cursor > 49){
                                        $text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result);
                                        $cursor = $cursor - 20;
                                } else {
                                        $cursor = 700;
                                        $page = $pdf->page();
                                        $page->mediabox('Letter');
                                        $font = $pdf->corefont('Helvetica-Bold');
                                        $text = $page->text();
                                        $text->font($font, 12);
                                        $text->translate(70, $cursor);
                                        $text->text($temp_result);
                                        $cursor = $cursor - 20;
                                } 
                        }
                
                        #printing closing line of query
                        $text = $page->text();
                        $text->font($font, 12);   
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $sth -> finish;
				
		} elsif ($fquery_bits[0] == 8){
			#querying database
                        my $sth = $dbh -> prepare("SELECT DISTINCT Cluster.motif_label FROM Cluster JOIN Word ON Word.cluster_name = Cluster.cluster_name AND Word.coverage = 1");
                        $sth -> execute();
			my @val = $sth->fetchrow_array();                                        

                        #printing header to PDF
                        my $cursor = 700;
                        $page = $pdf->page();
                        $page->mediabox('Letter');
                        $font = $pdf->corefont('Helvetica-Bold');
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * QUERY 8 * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $cursor = $cursor - 20;

			#printing query data to PDF
                        while(my @val = $sth->fetchrow_array()){

				if($cursor > 259){
                                        #setting cursor to correct location
                                        $cursor = $cursor - 200;
                                        my $temp_image = $val[0];
                                        my $temp_file = "temp.png";

                                        #saving image to a temporary .png, will be deleted later
                                        open(MYOUTFILE, ">$temp_file");
                                        print MYOUTFILE $temp_image;
                                        close(MYOUTFILE);
                                        
                                        #printing image to PDF
                                        my $gfx = $page -> gfx; 
                                        my $image = $pdf -> image_png($temp_file);
                                        $gfx -> image($image, 100, $cursor, 400, 200);
                                        $cursor = $cursor - 20;

                                        #deleting temporary .png
                                        unlink($temp_file);
                                } else {
                                        $cursor = 700;
                                        $page = $pdf->page();  

                                        #setting cursor to correct location
                                        $cursor = $cursor - 200;
                                        my $temp_image = $val[0];
                                        my $temp_file = "temp.png";

                                        #saving image to a temporary .png, will be deleted later
                                        open(MYOUTFILE, ">$temp_file");
                                        print MYOUTFILE $temp_image;
                                        close(MYOUTFILE);

                                        #printing image to PDF
                                        my $gfx = $page -> gfx;
                                        my $image = $pdf -> image_png($temp_file);
                                        $gfx -> image($image, 100, $cursor, 400, 200);
                                        $cursor = $cursor - 20;

                                        #deleting temporary .png 
                                        unlink($temp_file); 
                                }
                        }

                        #printing closing line of query
                        $text = $page->text();
                        $text->font($font, 12);
                        $text->translate(70, $cursor);
                        $text->text('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
                        $sth -> finish;
		} 
	}

#naming and saving report
$report_name = $report_name . ".pdf";
$pdf->saveas($report_name);
print "$report_name has been saved to your current directory.\n";
$dbh -> disconnect;

