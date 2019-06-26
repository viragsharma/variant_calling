#!/usr/bin/perl

## Virag Sharma, 2016. A variant calling pipeline that utilizes samtools, picard and GATK
## Based on "From FastQ data to high confidence variant calls: the GenomeAnalysis Toolkit best practices pipeline PMID: 25431634"
use strict;
use warnings;

my $usage   = "\nUSAGE: $0 bamFile ReferenceFasta (with a .fasta extension) outputFile(lists the variants in vcf format) knownIndels(vcf) dbSNP(vcf) \n\n";
my $purpose = "######\nPURPOSE: $0 is a variant calling pipeline which requires a bamFile and the referenceGenome fasta file\n\
$0 performs the following operations:
1. Sorts the input bam file
2. Removes duplicates
3. Creates dictionary and index for reference genome fasta file
4. Runs GATK-RealignTargetCreator on sites that require realignment and then realigns such sites using GATK-IndelRealigner
5. Performs base quality recalibration using information from knownIndels and dbSNP variants
6. Finally, runs GATK-HaplotypeCaller to report variants\n######\n";

die $usage.$purpose if (scalar(@ARGV) != 5);

my $bamFile     = $ARGV[0];
my $refFasta    = $ARGV[1];
my $outFile     = $ARGV[2];
my $knownIndels = $ARGV[3];
my $dbSNP       = $ARGV[4];

## Check if input parameters are correctly supplied
die "Input bamFile '$bamFile' does not exist\n" if (! -e $bamFile);
die "Reference genome fasta file '$refFasta' does not exist, or is not in the right format\n" if (! -e $refFasta || $refFasta !~/.fasta$/);
die "knownIndels file '$knownIndels' does not exist, or is not in the right format\n" if (! -e $knownIndels || $knownIndels !~/.vcf/);
die "dbSNP file '$dbSNP' does not exist, or is not in the right format\n" if (! -e $dbSNP || $dbSNP !~/.vcf/);
open(FO,">$outFile") || die "Error --> cannot write to the outFile '$outFile'\n";
close FO;

## Check if the relevant tools are installed and in the $PATH variable
checkTool("samtools");
checkTool("picard");
checkTool("gatk");

## All OK, time for business
## Get file name for the bamFile, the fasta file and the directory where the fasta file is present
my $fileName = `basename $bamFile`; chomp $fileName;
$fileName =~s/\.bam//;
my $dirNameFasta  = `dirname $refFasta`; chomp $dirNameFasta;
my $fileNameFasta = `basename $refFasta`; chomp $fileNameFasta;
$fileNameFasta =~s/.fasta//; 

## Step 0: Sort the bam file
my $sortedBam = "$fileName\_sorted";
my $sCall = "samtools sort $bamFile $sortedBam";
system($sCall) == 0 || die "Error sorting bam file, '$sCall' failed\n";
print "### Step1 done, bam file sorted ###\n\n";
$sortedBam = $sortedBam.".bam";

## Step 1: Remove duplicates with picard
my $sortedSam_NoDup = "$fileName\_NoDup.sam";
my $pCall = "picard MarkDuplicates INPUT=$sortedBam OUTPUT=$sortedSam_NoDup METRICS_FILE=metrics.txt";
system($pCall) == 0 || die "Error removing duplicates, '$pCall' failed\n";
print "### Step2 done, picardMarkDuplicates run successful ###\n\n";

my $refFastaDict = "$dirNameFasta/$fileNameFasta.dict";
## Step 2a: Create reference-fasta dictionary with CreateSequenceDictionary if it does not exist
if (! -e $refFastaDict){
	my $dCall = "picard CreateSequenceDictionary R=$refFasta O=$refFastaDict";
	system($dCall) == 0 || die "Error running CreateSequenceDictionary, '$dCall' failed\n";
	print "### Step2a done, picard CreateSequenceDictionary run successful ###\n\n";
}
## Step 2b: Create fasta index file
my $fastaIndex =  "$dirNameFasta/$fileNameFasta.fasta.fai";
if (! -e $fastaIndex){
	my $fiCall = "samtools faidx $dirNameFasta/$fileNameFasta.fasta";
	system($fiCall) == 0 || die "Error running samtools faidx, '$fiCall' failed\n";
	print "### Step2b done, samtools faidx run successful ###\n\n";
}
## Step 2c: Reorder the bamFile according to the reference file
my $sortedBamReordered = "$fileName\_NoDup_Reordered.bam";
my $reoCall = "picard ReorderSam I=$sortedSam_NoDup O=$sortedBamReordered REFERENCE=$refFasta";
system($reoCall) == 0 || die "Error running ReorderSam, '$reoCall' failed\n";
print " ## input bam file reordered according to the reference fasta file\n";

## And index this bam file
my $iCall1 = "samtools index -b $sortedBamReordered";
system($iCall1) == 0 || die "Error indexing bam file, '$iCall1' failed\n";

## Step 2: Run RealignTargetCreator
my $realignIntervals = "$outFile.intervals";
my $rCall1 = "gatk -T RealignerTargetCreator -R $refFasta -I $sortedBamReordered -known $knownIndels -o $realignIntervals";
system($rCall1) == 0 || die "Error running gatk-RealignerTargetCreator, '$rCall1' failed\n";
print "### Step2 done, RealignerTargetCreator run successful ###\n\n";

## Step 3: Now realign the stuff that RealignerTargetCreator has found
my $sortedBam_Realigned = "$fileName\_NoDup_Reordered_Realigned.bam";
my $rCall2 = "gatk -T IndelRealigner -R $refFasta -I $sortedBamReordered -targetIntervals $realignIntervals -known $knownIndels -o $sortedBam_Realigned";
system($rCall2) == 0 || die "Error running gatk-IndelRealigner, '$rCall2' failed\n";
print "### Step3 done, IndelRealigner run successful ###\n\n";
 
## Step 4: Performa base quality score recalibration
my $reCalDataFile  = "$fileName\_reCal.data";
my $recalibCall = "gatk -T BaseRecalibrator -R $refFasta -I $sortedBam_Realigned -knownSites $knownIndels -knownSites $dbSNP -o $reCalDataFile";
system($recalibCall) == 0 || die "Error running gatk-BaseRecalibrator, '$recalibCall' failed\n";
print "### Step4 done, BaseRecalibrator run successful ###\n\n";

## Step 5: Now apply the base-recalibration scores to the bam file
my $reCalBam = "$fileName\_NoDup_Reordered_Realigned_Recalibrated.bam";
my $printCall = "gatk -T PrintReads -R $refFasta -I $sortedBam_Realigned -BQSR $reCalDataFile -o $reCalBam";
system($printCall) == 0 || die "Error running gatk-PrintReads, '$printCall' failed\n";
print "### Step5 done, gatk-PrintReads run successful ###\n\n";
 
## Step 6: Create an index file
my $iCall2 = "samtools index -b $reCalBam";
system($iCall2) == 0 || die "Error indexing bam file, '$iCall2' failed\n";
print "### Step6 done, index created for the bam file '$reCalBam' ###\n\n";

## Step 7: Run variant caller:
my $vCall = "gatk -T HaplotypeCaller -R $refFasta -I $reCalBam -o $outFile";
system($vCall) == 0 || die "Error running gatk-HaplotypeCaller, '$vCall' failed\n";
print "\n\n\n#######\n$0 run successful, the input file was '$bamFile', the reference genome fasta file was '$refFasta' and the variants are written to '$outFile'\n";

sub checkTool {  ## This function checks if a tool that I need has been installed and is in '$PATH'
	my $tool = shift;
	my $status = "F";
	
	for my $path(split /:/,$ENV{PATH}) {
		if ( -f "$path/$tool" && -x _) {
			$status = "T";
			last;
		}
	}
	die "'$tool' not found, either it is not installed or is not in your \$PATH variable\n" if ($status eq "F");
}
