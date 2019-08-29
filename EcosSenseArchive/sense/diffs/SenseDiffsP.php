<?php
ini_set('memory_limit','2G'); 
// SCRIPT REQUIRES BASE AND DIFF FILES HAVE RESULTS IN THE SAME ORDER
//$BaseFile = "vulbase1.csv";
//$DiffFile = "vulFfishery1.csv";
$DiffFile = $argv[1];
$BaseFile = $argv[2];

//echo "$argv[1]\n"; exit(0);
  $Bbase = array();
  $DiffSeries = array();
  $TestS = -1;
  $Scenes = 0;
  $handle   = fopen($BaseFile,'r');
  $handle2  = fopen($DiffFile,'r');
      $firstline = fgetcsv($handle,5000);  
      $f2        = fgetcsv($handle2,5000); 
			while ($line=fgetcsv($handle,5000)){
			       $l2  =fgetcsv($handle2,5000);
            for ($i=0; $i<count($firstline); $i++){$$firstline[$i] = $line[$i];}
            for ($i=0; $i<count($firstline); $i++){$fn="$firstline[$i]_2"; 
                                                   $$fn = $l2[$i];}
            $key = "$Seed,$System,$Species,$Year";
            $k2  = "$Seed_2,$System_2,$Species_2,$Year_2";
            if ($key != $k2){echo "$key : $k2 Oops out of order!\n"; exit(1);}
            //if ($outBB<=0.0){echo "$key Base biomass 0!\n"; exit(1);}
						//$Bbase[$key] = $outBB;
            if (($outBB*$outPB)>0.0){$diff = (($outBB_2*$outPB_2)-($outBB*$outPB))/($outBB*$outPB);}
            else           {$diff = 0.0;}
            //echo "$System $outBB $outBB_2 $diff\n";
            $varKey = "$Species,$Year";
            if (!array_key_exists($varKey,$DiffSeries)){$DiffSeries[$varKey]=array();}
            $DiffSeries[$varKey][] = $diff;  // This stores is the next slot in an array
            if ($System != $TestS){$Scenes++; fprintf(STDERR, "Run $Scenes (Seed $Seed : Systen $System)\n");}
            $TestS = $System;
            
			 }
       fclose($handle);
       fclose($handle2);
       
       echo "BaseScenario,DiffScenario,Species,Year,Nsystems,";
       for ($i=0; $i<=100; $i+=5){echo "CI_$i,";} echo "\n";
       
       foreach ($DiffSeries as $DiffKey => $DiffSet){
			    sort($DiffSet);
			    $NN = count($DiffSet)-1;
			    $NS = count($DiffSet);
					echo "$BaseFile,$DiffFile,$DiffKey,$NS,";
					for ($i=0; $i<=100; $i+=5){
					     $ind = (int)(0.01*$i*$NN);
					     echo "$DiffSet[$ind],";
					}
					echo "\n";
			 }
       
?>
