package temp;

import java.util.ArrayList;

/*******************************************************
 * Copyright (C) 2018 Pratyusha Pogaru <plpratyusha@gmail.com>
 * 
 * This file is part of the CS 123B final project for Spring '18.
 * 
 * TrimAdvisor cannot be copied and/or distributed without the express
 * permission of Pratyusha Pogaru.
 *******************************************************/

public class TrimAdvisor {

	public static void main(String[] args) {
		
		String fileToParse = ">sp|P53347|ONCM_MOUSE Oncostatin-M OS=Mus musculus OX=10090 GN=Osm PE=1 SV=1\n" + 
				"------------------MQTRLLRTLLSLTLSLLILSMA--LANRGCSNSSSQLLSQLQNQANLTGNTESLLEPYI\n" + 
				"RLQNLNTPDLRAACTQHSVAFPSEDTLRQLSKPHFLSTVYTTLDRVLYQLDALRQKFLKT\n" + 
				"PAFP--KLDSARHNILGIRNNVFCMARLLNHSLEIPEPTQTDSGASR--S-T\n" + 
				"TTPDVFNTKIGSCGFLWGYHRFMGSVGRVFREWDDGSTRSRRQSPLRARRKGTRRIRVRH\n" + 
				"KGTRRIRVRRKGTRRIWVRRKGSRKIRPSRSTQSPTTRA\n" + 
				">sp|Q65Z15|ONCM_RAT Oncostatin-M OS=Rattus norvegicus OX=10116 GN=Osm PE=2 SV=1\n" + 
				"---------------MRAQPPPRTLLSLALALLFLSMS--WAKRGCSSSSPKLLSQLKSQANITGNTASLLEPYI\n" + 
				"LHQNLNTLTLRAACTEHPVAFPSEDMLRQLSKPDFLSTVHATLGRVWHQLGAFRQQFPKI\n" + 
				"QDFP----ELERARQNIQGIRNNVYCMARLLHPPLEIPEPTQADSGTSR--PTT\n" + 
				"TAPGIFQIKIDSCRFLWGYHRFMGSVGRVFEEWGDGSRRSRRHSPLWAWLKGDHRIRPSR\n" + 
				"SSQSAMLRS-LVPR------\n" + 
				">tr|A0A0D9RAW5|A0A0D9RAW5_CHLSB Oncostatin M OS=Chlorocebus sabaeus OX=60711 GN=OSM PE=4 SV=1\n" + 
				"--------------------------MGVPLTRRTLLSLVLALLFPSMASMAAMGSCSKEYRMLLGQLQKQTDLMQDTSRLLDPYI\n" + 
				"RIQGLDIPKLREHCRESPGAFPSEETLRGLGRRGFLQTLNDTLGRVLHRLADLEQHLPKA\n" + 
				"QDLERSGLNIEDLEKLQMARPNVLGLRNNIYCMAQLLDN-SDMTEPTKAGRGAPQPPTPT\n" + 
				"PTSDVFQRKLEGCSFLHGYHRFMHSVGQVFSKWGESPNRSRRHSPHQALRKGMRRTRPSR\n" + 
				"KGNRLMPR-GQLPR-----\n" + 
				">tr|E2RQT2|E2RQT2_CANLF Oncostatin M OS=Canis lupus familiaris OX=9615 GN=OSM PE=4 SV=1\n" + 
				"---------------------MQAQLLWRTLPSLVLGLLFLSM---PAMGSCSDKYPELLGQLQKQADFMQHTNTLLDLYI\n" + 
				"RSQGLDKNGLKEHCRERPGAFPSKDALQRLSRRVFLRTLDTTLGQVLLRLAALEQDIPKA\n" + 
				"QDLE----------MLSGVKLNIRGFKNNIHCMAQLLPGSSETTEPTPTSPGASP--SPT\n" + 
				"PTLDTFQRRLEGCRFLHGYHRFMRSVGQVFREWGKSLSRSRRHSPHQGLLKGARRMQLSG\n" + 
				"RNKRLMPR-VQLAPGSHRGAPGG-----\n" + 
				">tr|M3YP54|M3YP54_MUSPF Oncostatin M OS=Mustela putorius furo OX=9669 GN=OSM PE=4 SV=1\n" + 
				"-------------------------RGRCMRSRTVLGLVLGLLFLST---PVMGSCLDKDQELLRQLQKQADIMQRTSMLLNPYI\n" + 
				"QSQGLDKDGLKEHCRERPGTFPSKEALQRLSRQELLRALNTMLDHVLHRLMALQQDIPKA\n" + 
				"QDLE----------TVNRAKQNIRGFKNNIHCMAQLLPGSSERTEPPPIGPGTSP--SPT\n" + 
				"PTPDTFQRRLEGCRFLHGYHRFMHSVGQVFREWAESPSRSRRHSPRRGLWKGARRVPLSR\n" + 
				"GNKRLLLR-GQLPR----------"; //Directly paste sequence from an aligned training set into the quotes. (Place cursor between quotes and paste)

		//Creates variables and data structures
		int fileLength = fileToParse.length();
		int numFastas = 0;
		int startFasta = 0; //index of starting > in a fasta
		int endFasta = 0; //index of last char in a fasta
		ArrayList<Integer> startF = new ArrayList<Integer>();
		ArrayList<Integer> endF = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> fastaBounds = new ArrayList<ArrayList<Integer>>();
		ArrayList<String> fastaFiles = new ArrayList<String>();
		ArrayList<Integer> indelCount = new ArrayList<Integer>();
		ArrayList<Integer> longestIndelRun = new ArrayList<Integer>();
		ArrayList<Integer> startIndel = new ArrayList<Integer>();
		ArrayList<Integer> endIndel = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> longestIndelBounds = new ArrayList<ArrayList<Integer>>();

		
		//Prints error message if input file is empty
		if (fileLength == 0) {
			System.out.println("file empty");
		}

		
		//Calculates start and end indices of each fasta file; also counts # fasta files total
		for (int n = 0; n < fileLength; n++) {
			if (fileToParse.substring(n, n+1).equals(">")) {				

				int oldStartFasta = startFasta;
				startFasta = n;
				startF.add(numFastas, startFasta);

				if (numFastas > 0) {
					int length = startFasta - oldStartFasta - 1;
					endFasta = oldStartFasta + length;
					endF.add(numFastas-1, endFasta);
				}
				numFastas++;
			}
		}
		endF.add((fileLength-1));
		System.out.println("Total number of fasta files: " + numFastas);
		System.out.println();
		

		//Gets start/end indices of each fasta and adds to fastaBounds ArrayList
		for (int i = 0; i < startF.size(); i++) {
			ArrayList<Integer> temp = new ArrayList<Integer>();
			temp.add(startF.get(i));
			temp.add(endF.get(i) + 1); //changed to + 1 to account for last chars

			fastaBounds.add(i, temp);
		}

		
		//Gets each fasta file separately and adds to fastaFiles ArrayList
		for (int i = 0; i < numFastas; i++) {
			String temp = fileToParse.substring(fastaBounds.get(i).get(0), fastaBounds.get(i).get(1));
			fastaFiles.add(temp);
		}


		//Deletes first line of each fasta
		for (int i = 0; i < fastaFiles.size(); i++) {
			String f = fastaFiles.get(i);
			fastaFiles.set(i, f.substring(f.indexOf('\n')+1));
		}
		
		//Counts indels of each fasta and populates indelCount ArrayList
		for (int i = 0; i < fastaFiles.size(); i++) {
			int count = 0;
			String currentFasta = fastaFiles.get(i);

			for (int j = 0; j < currentFasta.length(); j++) {
				if (currentFasta.substring(j, j+1).equals("-")) {
					count++;;
				}
			}
			indelCount.add(count);
		}

		
		//Also counts longest run of indels for each fasta and populates longestIndelRun ArrayList
		for (int i = 0; i < fastaFiles.size(); i++) {
			int maxRun = 0; 
			int theMax = 0;
			int startIndex = 0;
			int finalStartIndex = 0;
			int endIndex = 0;
			String currentFasta = fastaFiles.get(i);
			
			for (int j = 0; j < currentFasta.length(); j++) {
				if (j+2 <= currentFasta.length() && (currentFasta.substring(j, j+1).equals("-")) && (currentFasta.substring(j+1, j+2).equals("-"))) {
					maxRun++;
					startIndex = j - maxRun;
				}
				else {
					if (maxRun + 1 > theMax) {
						theMax = maxRun + 1;
						finalStartIndex = startIndex + 1;
						endIndex = finalStartIndex + theMax;
					}
					maxRun = 0;
					startIndex = 0;
				}
				
			}
			longestIndelRun.add(theMax);
			startIndel.add(finalStartIndex);
			endIndel.add(endIndex);
		}
		
		
		//Gets indices of the bounds of the longest indel run for each fasta. Populates longestIndelBounds ArrayList
		for (int i = 0; i < startIndel.size(); i++) {
			ArrayList<Integer> temp = new ArrayList<Integer>();
			temp.add(startIndel.get(i));
			temp.add(endIndel.get(i));
			
			longestIndelBounds.add(i, temp);
		}
		
		
		//Prints each fasta sequence & corresponding information
		for (int i = 0; i < fastaFiles.size(); i++) {
			System.out.println("Sequence for fasta #" + i); 
			System.out.println(fastaFiles.get(i));
			System.out.println("Number of indels in this fasta: " + indelCount.get(i));
			System.out.println("Maximum run of indels in this fasta: " + longestIndelRun.get(i));
			System.out.println("Longest indel bounds: (" + longestIndelBounds.get(i).get(0) + ", " + longestIndelBounds.get(i).get(1) + ")");
			System.out.println("\n");
		}

		
		//Calculates if trimming is necessary
		int countEae = 0; //number of fastas that End At the End of the sequence
		int countSas = 0; //number of fastas that Start At the Start of the sequence
		
		for (int i = 0; i < fastaFiles.size(); i++) {
			if (longestIndelBounds.get(i).get(0) == 0) {
				countSas++;
			}
			if (longestIndelBounds.get(i).get(1) == fastaFiles.get(i).length()) {
				countEae++;
			}
		}
		
		
		//Prints message advising trimming if needed
		if ( (countSas >= fastaFiles.size()/4 && countSas < fastaFiles.size()/4) || (countEae >= fastaFiles.size()/4 && countEae < fastaFiles.size()/4) ) {
			System.out.println("Based on the positioning of the indels in your aligned fasta sequences, it might be a good idea to trim your sequences.");
		}
		else if (countSas >= fastaFiles.size()/2) {
			System.out.println("Many sequences (half or more) in your alignment file start with indels. Trimming your sequences down is advised.");
		}
		else if (countEae >= fastaFiles.size()/2) {
			System.out.println("Many sequences (half or more) in your alignment file end with indels. Trimming your sequences down is advised.");
		}
		
		System.out.println();
		System.out.println("---FINAL VERSION---");
		

	}

}
