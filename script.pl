#!/usr/bin/perl
use strict;
use CGI;

my $start_run = time();
my $date = `date`."\n";

$CGI::POST_MAX = 1024 * 3000000;			#5Mb maximum upload
my $random = int(rand(1000000000000000));	#Random number from 0 to 999999999999999
my $formulario = new CGI;					#GET INITIAL INFO FROM INDEX.HTML##
print $formulario->header;

#######################
#CREATE FOLDERS & FILES
#######################
my $upload_dir = "/var/www/html/ab/results";
#my $db = "/home/databases/ab/cel.fasta";
#my $db= "/home/databases/anablast/Uniref50/04-2014/uniref50.fasta"; ###para DROME
#my $db=  "/home/databases/anablast/Uniref50/01-2016/uniRef50.fasta;  ###para DROME
#my $db=  "/home/databases/anablast/Uniref50/uniref50_mod/Uniref50_mod.fasta";
my $db=  "/home/databases/anablast/Uniref50/05-2017/uniref50.fasta"; ###para DROME
#my $db=  "/home/databases/anablast/Uniref50/uniref50_mod/uniref50.fasta"; ###original
#my $db = "/home/databases/anablast/Uniref50/01-2012/uniref50.fasta"; ###para DROME
my $results_web = "http://www.bioinfocabd.upo.es/ab/results/$random";	
my $random_dir = "$upload_dir/$random";
`mkdir $random_dir`;
`chmod 777 $random_dir`;	
my $output1 = "$random_dir/$random-output1.blast";

&logfile ("$random_dir/script.log", $date);

#close STDIN; close STDERR; close STDOUT;

###############
#GET INPUT DATA
###############
my $sequences = $formulario->param("file1");
my $blou = $formulario->param("file2");																
my $pasted = $formulario->param("seqs_copy");
my $Bscore= $formulario->param("Bit_Score");
if ((!$sequences) && (!$pasted) && (!$blou)) {
	my $error = "To start the execution of AnABlast, please use as input one sequence.";
	&error($error);
}
my $file = "$random_dir/$random.fasta";

#UPLOADED FILE
if ($sequences) {		
	my $upload_filehandle = $formulario->upload("file1");	
	open UPLOADFILE, ">$file" or die $!;
	my $f = 0;
	while (<$upload_filehandle>) {
		if ($f == 0) {
			$f = 1;
			if ($_ !~ /^>/) {
				print UPLOADFILE ">$random\n"; #If input is in a file, with or without headline
			}
			elsif ($_ =~ /^>/) {
				$_ =~ s/(>.+?)[^a-zA-Z0-9].+/$1/;
			}
		}
		print UPLOADFILE $_;
	}
	close UPLOADFILE;
}		
#PASTED INPUT
if ($pasted) {			
	if ($pasted !~ /^>/) {
		$pasted = ">$random\n".$pasted; #If input is pasted, with or without headline
	}
	else {
		$pasted =~ s/(>.+?)[^a-zA-Z0-9].+/$1/m; # chenge seq id due to extrange characters
	}
	open COPYFILE, ">$file" or die $!;
	print COPYFILE "$pasted\n";
	close COPYFILE;
}		

#UPLOADED BLAST FILE
if ($blou) {
	my $upload_filehandle = $formulario->upload("file2");	
	open UPLOADFILE, ">$output1" or die $!;
	while (<$upload_filehandle>) {
		print UPLOADFILE $_;
	}
	close UPLOADFILE;	
}

###############################
#PARSE FASTA FILE AND CHECK IT#
###############################
my $fasta_id;
my $fasta_length;
my $fasta_nn;
open FASTA, $file or die $!;
while (my $line = <FASTA>) {
	$line =~ s/[\n\r]+$//g;
	next if ($line =~ m/^$/);
	if ($line =~ m/^>([^\s\t\n]+)/) {
		$fasta_id = $1;
		next;
	}
	else {
		$fasta_length += length($line);
		$fasta_nn .= $line;
	}
}
close FASTA or die $!;

if ($blou) {
	if ($fasta_length > 1000000) { #Seq longer than 1Mb
		my $error = "To start the execution of AnABlast, please use as input a sequence shorter than 1 Mb.";
		&error($error);	
	}
}
else {
	if ($fasta_length > 50000) { #Seq longer than 25000
		my $error = "To start the execution of AnABlast, please use as input a sequence shorter than 50KB.";
		&error($error);	
	}
}

if ($fasta_nn =~ m/[^ACGTUNacgtun]/) { #Seq not in DNA alphabet
	my $error = "To start the execution of AnABlast, please use as input a nucleotide sequence!";
	&error($error);
}

&logfile ("$random_dir/script.log", "sequence id: $fasta_id\nsequence length: $fasta_length\n\n");

###########
#RUN QUEUE#
###########
if (fork) {
	&summary($random,$results_web);
	&temporary($random,$random_dir,$results_web, "Please, be patient, job has been queued up.");
	exit;
}

my $f = 0;
while ($f == 0) {
	my @ps = `ps -eF | grep -iE "blastx|sma3s" | grep -v grep`;
	my $n = @ps; 
		if ($n == 0) {
			$f = 1;
			last;
		}
	sleep 10;
}

if (fork) {
	&temporary($random,$random_dir,$results_web, "Please, be patient, it will take a short while.");
	exit;
}

#TBLASTN
#`makeblastdb -in $file -dbtype nucl`;
#`tblastn -db $file -query $db -evalue 100 -outfmt '6 qseqid sseqid sstart send evalue bitscore sframe' -matrix BLOSUM90 -seg yes -max_target_seqs 10000 -out $output1 -num_threads 8` if (!$blou); 

#BLASTX
my $evalue = 200;
my $max_target_seqs = 10000;
my $seg = "no";
&logfile ("$random_dir/script.log", "\n\nRUNNING BLAST:\nBlast parameters: evalue = $evalue\t max_target_seqs= $max_target_seqs\n\n");

`blastx -db $db -query $file -evalue $evalue -outfmt '6 sseqid qseqid qstart qend evalue bitscore qframe' -matrix BLOSUM90 -seg $seg -max_target_seqs $max_target_seqs -out $output1 -num_threads 8 &>> $random_dir/script.log` if (!$blou);

#BLASTX ORIGINAL
#`blastx -db $db -query $file -evalue 11500 -outfmt '6 sseqid qseqid qstart qend evalue bitscore qframe' -matrix BLOSUM90 -seg yes -max_target_seqs 10000 -out $output1 -num_threads 8 &> $random_dir/blastx.log` if (!$blou);


#BLASTX
#`blastx -db $db -query $file -evalue 50000 -outfmt '6 sseqid qseqid qstart qend evalue bitscore qframe' -matrix BLOSUM90 -seg no -max_target_seqs 100000000 -out $output1 -num_threads 8` if (!$blou);


if (-z $output1) { #blast report is empty
	my $error = "Blast report is empty, please check your input!";
	&error2($error, $random_dir, $random);
}
else { # check blast and sequence id
	open BLAST, $output1 or die $!; 
	while (<BLAST>) {
		if ($_ !~ m/^.+?\t.+?\t\d+\t\d+\t.+?\t.+?\t-?\d$/) {
			my $error = "Blast report has wrong format, please check your input!";
			&error2($error, $random_dir, $random);
			
		} 
		my ($c1, $c2) = split(/\t/, $_);
		if ($c2 ne $fasta_id) {
			my $error = "Blast report has wrong sequence id, please check your input!";
			&error2($error, $random_dir, $random);
		}
		last;
	}
	close BLAST or die $!;
}
#######################
#AB INITIO PREDICTIONS#
#######################
&logfile ("$random_dir/script.log", "\n\nRUNNING AB INITIO PREDICTIONS:\n");

my $n = 0;
my $id;
open BLAST, $output1 or die $!;

while (<BLAST>) {
	my @colum = split(/\t|_/, $_);
	if (length $colum[1] == 6 and $colum[6] > $n) {
		$n = $colum[6];
		$id = $colum[1];
	}
}
close BLAST or die $!;
open TAX, ">$random_dir/best_tax.txt" or die $!;
print TAX $id;
close TAX or die $!;

`wget http://www.uniprot.org/uniprot/$id.txt -O $random_dir/taxonomy.dat &>> $random_dir/script.log`;

open AUG, "./augustus_taxonomy.tsv" or die $!;  
my @augustus_tax = <AUG>;
close AUG or die $!;	

my $species = "human"; # por defecto se ejecuta augustus con human como parametro
my $oc;
my $os;
my $a = 0; 
my $p = 0;
open TAX, "$random_dir/taxonomy.dat" or die $!;
while (<TAX>) {
	chomp;
	if ($_ =~ m/^OC\s{3}(.+)/) {
		$oc .= $1;
	}
	elsif ($_ =~ m/^OS\s{3}(.+)/) {
		$os .= $1;
	}
}
close TAX or die $!;

if ($oc =~ m/bacteria/i) {
	$a = 0;
	$p = 1;
}
else {
	foreach my $line (@augustus_tax) {
		next if ($line =~ m/^#/);
		my ($identifier, $sp, $tax) = split(/\t/, $line);
		if ($os =~ m/$sp/i) {
			$species = $identifier;
			$a = 1;
			$p = 0;
			last;
		}
	}
}

if ($a == 0 and $p == 0) {
	my @oc = split(/;/, $oc);
	foreach (@oc) {
		chomp;
		$_ =~ s/^\s+|\s+$//g;
		foreach my $line (@augustus_tax) {
			next if ($line =~ m/^#/);
			my ($identifier, $sp, $tax) = split(/\t/, $line);
			next if (! $tax);
			my @taxon = split(/;/, $tax); 
			for (my $i = 0; $i <= $#taxon; $i++) {
				chomp $taxon[$i];
				$taxon[$i] =~ s/^\s+|\s+$//g;
				if ($_ =~ m/$taxon[$i]/) {
					$species = $identifier;
					last;
				}			
			}
		}	
	}
}
if ($a == 0 and $p == 0) {
	$a = 1; # si no encuenta ningún organismo cercano se ejecuta augustus como human
}

`/home/apps/augustus.2.5.5/bin/augustus --species=$species --gff3=on --AUGUSTUS_CONFIG_PATH=/home/apps/augustus.2.5.5/config/ $file > $random_dir/augustus.gff` if ($a == 1);
`prodigal -p meta -i $file -f gff -o $random_dir/prodigal.gff &>> $random_dir/script.log` if ($p == 1);

##########
#ANABLAST#
##########
&logfile ("$random_dir/script.log", "\n\nRUNNING ANABLAST SCRIPT:\n");
`./anablast_v1.1.pl $file $output1 $random_dir . $Bscore . &>> $random_dir/script.log`;

#########
#RESULTS#
#########

#Build jbrowse dataset
my $jb_dir = "/var/www/html/jbrowse_anablast/JBrowse-1.12.1"; 
&add_to_jbrowse ($random, $random_dir, $jb_dir);
#Build result web
&result($random, $random_dir);	
`cp $random_dir/tmp.html $random_dir/results_$random.html`;

$date = `date`;
my $end_run = time();
my $run_time = $end_run - $start_run;
my $sec_kb = (1000 * $run_time) / $fasta_length;
&logfile ("$random_dir/script.log", "\n\n########### SUMMARY ###########\n$date\nRun time: $run_time seconds\nSequence length: $fasta_length\nSeconds/kb: $sec_kb\nBlast evalue: $evalue\nBlast max_target_seqs: $max_target_seqs\nBlast seg: $seg\n");

exit;
					
#############
#SUBROUTINES#
#############
sub logfile () {
	my ($file, $message) = @_;
	open LOG, ">>$file" or die $!;
	print LOG $message;
	close LOG or die $!;
}

sub result () {

	my $rand = shift;
	my $rdir = shift;
	my $html = "http://www.bioinfocabd.upo.es/jbrowse_anablast/JBrowse-1.12.1/index.html?data=datasets%2F".$rand."&amp;tracks=DNA%2CMultiBigWig%2Cpeaks%2CORF%2CProdigal%2CAugustus&loc=$fasta_id%3A1..$fasta_length";

	open TMP, ">$rdir/tmp.html" or die $!;

	print TMP <<EOF;
 <html lang="en">
  <head>
    <meta content="text/html; charset=UTF-8" http-equiv="content-type">
    <link rel="shortcut icon" href="/ab/favicon.ico">
    <title>AnABlast</title>
    <script>function abrir(url) { open(url,'','top=100,left=700,width=500,height=500') ; }</script>
    <style type="text/css" media="screen">
	.tab {		padding-left:5em								}
	.bod {		font-family:arial;							}
	form {    display: inline-block;   margin-left: auto;   margin-right: auto;   text-align: left; }
	</style>
  </head>
  <body class="bod"> 
	<table style="border: 1px solid; width: 100%;" align="center" border="0" cellpadding="0" cellspacing="0" align="center">
      <tbody> 
        <tr> 
          <td style="text-align: left; width: 55px;"> <a href="/ab/"><img longdesc="/ab/logo_ab.png" 
              title="logo" alt="logo" src="/ab/logo_ab.png" height="50" width="50"></a><br> 
          </td> 
          <td style="text-align: left; width: 150px;"><br> 
            <h2 style="heigth: 12px; margin-top: -6px; margin-bottom : 12px;">AnABlast</h2> 
          </td> 
          <td style="font-size:14px">  
          <u>Result tips:</u> <i><font color="#202020">Select tracks from the left block to activate views <font color=red>|</font> 
             Push on the items to expand information <font color=red>|</font> 
             Open drop-down menu to download sequences and more advanced settings</font></i>
          </td> 
          <td style="text-align: center; width: 20px;"> 
            <h3 style=" text-align: right; margin-right: 10px; margin-top: 12px;"> <a target="_blank" href="http://www.bioinfocabd.upo.es/ab/help.html">Help</a></h3>  
          </td> 
        </tr> 
      </tbody>
    </table>
    <div style="margin: 0 auto;"> <iframe style="border: 1px solid black" src="$html"
        height="100%" width="100%"> </iframe> </div>
    <br>
  </body>
 </html>
EOF

close TMP or die $!;
}

sub error () {
	my ($error) = @_;
	&logfile ("$random_dir/script.log", "\n$error\n");
	print <<EOF;
	<table width=\"900\" cellspacing=\"0\" cellpadding=\"10\" border=\"0\" align="center">
	<tr><td style=\"border-style: dotted\" align=\"center\" color="#FFFFFF" size="5"><br>$error<br><br>
	</td></tr></table></body><p>
EOF
exit;
}

sub error2 () {
	my ($error, $random_dir, $random) = @_;
	&logfile ("$random_dir/script.log", "\n$error\n");
	open (TMP, ">$random_dir/results_$random.html") or die $!;
	print TMP <<EOF;
	<table width=\"900\" cellspacing=\"0\" cellpadding=\"10\" border=\"0\" align="center">
	<tr><td style=\"border-style: dotted\" align=\"center\" color="#FFFFFF" size="5"><br>$error<br><br>
	</td></tr></table></body><p>
EOF
	close TMP or die $!;
	exit;
}

sub add_to_jbrowse {
	#
	#
	my $random = shift;
	my $random_dir = shift;
	my $jb_dir = shift;
	my $jb_random_dir = "$jb_dir/datasets/$random";
	$jb_random_dir =~ m/(\/jbrowse.*)$/;
	my $dataset_dir = $1;

	#create random directory
	`mkdir $jb_random_dir`;
	#load reference sequence into jbrowse
	`$jb_dir/bin/prepare-refseqs.pl --fasta $random_dir/$random.fasta --out $jb_random_dir --key "Nucleotide sequence" --seqType dna`;
	#create new dataset
	open CONF, ">>$jb_dir/jbrowse.conf" or die $!;
	print CONF <<DATASET;

[datasets.$random]
url  = ?data=datasets/$random
name = $random

DATASET
	close CONF or die $!;

	`cp $random_dir/*.bw $jb_random_dir`;

	open CONF, ">>$jb_random_dir/tracks.conf" or die $!;
	print CONF <<MWG;

[tracks.MultiBigWig]
key=AnABlast
type=MultiBigWig/View/Track/MultiWiggle/MultiXYPlot
style.height=200
storeClass=MultiBigWig/Store/SeqFeature/MultiBigWig
max_score=500
autoscale=global
category=AnABlast
description=AnABlast: Non-significan aligment accumulation. In differents shades of green accumulations over fordware frames strand an in differents shades of red accumulations over reverse frames strand"
MWG
	
	if (-e "$jb_random_dir/frame1.bw") {
		print CONF <<MWG
urlTemplates+=json:{"url":"$dataset_dir/frame1.bw", "name": "frame1", "nonCont" : true, "fill": "true", "color": "rgba(0, 230, 0, 0.4)"}
MWG
	}
	if (-e "$jb_random_dir/frame2.bw") {
		print CONF <<MWG
urlTemplates+=json:{"url":"$dataset_dir/frame2.bw", "name": "frame2", "nonCont" : true, "fill": "true", "color": "rgba(0, 177, 0, 0.4)"}
MWG
	}
	if (-e "$jb_random_dir/frame3.bw") {
		print CONF <<MWG
urlTemplates+=json:{"url":"$dataset_dir/frame3.bw", "name": "frame3", "nonCont" : true, "fill": "true", "color": "rgba(0, 100, 0, 0.4)"}
MWG
	}
	if (-e "$jb_random_dir/frame-1.bw") {
		print CONF <<MWG
urlTemplates+=json:{"url":"$dataset_dir/frame-1.bw", "name": "frame-1", "nonCont" : true, "fill": "true", "color": "rgba(255, 0, 0, 0.4)"}
MWG
	}
	if (-e "$jb_random_dir/frame-2.bw") {
		print CONF <<MWG
urlTemplates+=json:{"url":"$dataset_dir/frame-2.bw", "name": "frame-2", "nonCont" : true, "fill": "true", "color": "rgba(177, 0, 0, 0.4)"}
MWG
	}
	if (-e "$jb_random_dir/frame-3.bw") {
		print CONF <<MWG
urlTemplates+=json:{"url":"$dataset_dir/frame-3.bw", "name": "frame-3", "nonCont" : true, "fill": "true", "color": "rgba(100, 0, 0, 0.4)"}
MWG
	}
	
	# add peaks and ORF tracks
	`$jb_dir/bin/flatfile-to-json.pl --gff $random_dir/peaks.gff --out $jb_random_dir --trackLabel peaks --key Peaks --Tracktype JBrowse/View/Track/CanvasFeatures --clientConfig '{"color" : "violet"}' --metadata '{"description": "Peaks: Click on item to show more information as annot and top blast", "category": "AnABlast" }'`;
	`$jb_dir/bin/flatfile-to-json.pl --gff $random_dir/orf.gff --out $jb_random_dir --trackLabel ORF --key "Open Reading Frames" --Tracktype JBrowse/View/Track/CanvasFeatures --clientConfig '{"color" : "blue"}' --metadata '{"description": "ORF: Open Reading frames (start codon to stop codon) inside Peaks frame", "category": "Gene Finders" }'`;
	`sed -E "s/\ttranscript\t/\tmRNA\t/" $random_dir/augustus.gff > $random_dir/augustus2.gff` if (-e "$random_dir/augustus.gff");
	`$jb_dir/bin/flatfile-to-json.pl --gff $random_dir/augustus2.gff --out $jb_random_dir --trackLabel Augustus --key 'Eukaryotic prediction (Augustus)' --Tracktype JBrowse/View/Track/CanvasFeatures --clientConfig '{"color" : "orange"}' --metadata '{"description": "Augustus: ab initio eukaryotic gene prediction with Augustus program", "category": "Gene Finders" }'` if (-e "$random_dir/augustus.gff");
	`$jb_dir/bin/flatfile-to-json.pl --gff $random_dir/prodigal.gff --out $jb_random_dir --trackLabel Prodigal --key 'Prokaryotic prediction (Prodigal)' --Tracktype JBrowse/View/Track/CanvasFeatures --clientConfig '{"color" : "red"}' --metadata '{"description": "Prodigal: ab initio prokaryotic gene prediction with Prodigal program", "category": "Gene Finders" }'`;
	`$jb_dir/bin/json2conf.pl $jb_random_dir/trackList.json >> $jb_random_dir/tracks.conf`;
	`echo "{}" > $jb_random_dir/trackList.json`;

	close CONF or die $!;
}

sub summary () {
  my ($random,$web) = @_;
  print "<html><head><title>Results in short time</title>\n";
  print "<meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1'>\n";
  print "<META HTTP-EQUIV='Pragma' CONTENT='no-cache'>\n";
  print "<SCRIPT>setTimeout(\"location='$web/results_$random.html?cgi='+Math.random()\",1000);</SCRIPT></head></html>\n";
}

sub temporary () {
  my ($random,$random_dir,$web, $message) = @_;
  open (TMP, ">$random_dir/results_$random.html") or die $!;
  print TMP "<html><head><title>Results in short time</title>\n";
  print TMP "<meta http-equiv='Content-Type' content='text/html; charset=iso-8859-1'>\n";
  print TMP "<META HTTP-EQUIV='Pragma' CONTENT='no-cache'><script>\n";
  print TMP "setTimeout(\"location='$web/results_$random.html?cgi='+Math.random()\",10000);</script></head>\n";
  print TMP "<body><p>&nbsp;</p><p>&nbsp;\n";
  print TMP "<div align='center' style='border-width: 1px; border-style: dashed; '><b>\n";
  print TMP "<p>This page will be reloaded in 10 seconds.</p>\n";
  print TMP "<p><img src='/images/running.gif' alt='running' /></p>\n";
  print TMP "<p>$message<p><script>\n";
  print TMP "document.write('<a href=$web/results_$random.html'+'?cgi='+Math.random()+'>$web/results_$random.html</a>');\n";
  print TMP "</script></b></div>";
  #print TMP "</body></html>\n";
  close TMP or die $!;;
}
