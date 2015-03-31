
/**
 *
 * @author Kanishk Asthana kasthana@eng.ucsd.edu GitHub: github.com/kanishkasthana
 * Created on March 25, 2015 to filter ChipSeq Peaks generated by Homer
 */

import java.io.*;
import java.util.*;

public class homerPeakFilter {

    /**
     * @param args the command line arguments
     */
    
    public static List merge(List <Double>list1,List <Double>list2){
        List <Double>sortedList=new <Double>ArrayList();
        int firstPos=0,secondPos=0;
        Double current1,current2;
        
        while(firstPos<list1.size() && secondPos<list2.size()){
            current1=list1.get(firstPos);
            current2=list2.get(secondPos);
            
            if(current1<current2){
                sortedList.add(list1.get(firstPos));
                firstPos++;
            }
            
            else{
                sortedList.add(list2.get(secondPos));
                secondPos++;
            }
        }
        
        for(int i=firstPos;i<list1.size();i++){
            sortedList.add(list1.get(i));
        }
        
        
        for(int i=secondPos;i<list2.size();i++){
            sortedList.add(list2.get(i));
        }
        
        return sortedList;
    }
    
    public static List mergeSort(List<Double> list){
        
        if(list.size()==1){
            return list;
        }
        int firstHalfIndex=list.size()/2;
        List <Double>firstHalf=list.subList(0, firstHalfIndex);
        List <Double>secondHalf=list.subList(firstHalfIndex,list.size());
        List <Double>sortedFirstHalf=mergeSort(firstHalf);
        List <Double>sortedSecondHalf=mergeSort(secondHalf);
        List <Double>sortedList=merge(sortedFirstHalf,sortedSecondHalf);
        
        return sortedList;
    }
        
    public static void main(String[] args) {
        
        String inputFileName=args[0];
        String outputFileName=args[1];
        
        //Stored the peaks above this percentile
        int filterLevel=Integer.parseInt(args[2]);
        
        try{
            
            PrintWriter out= new PrintWriter(new FileWriter(outputFileName));
            List<String> inputs= new <String>ArrayList();
            File newFile=new File(inputFileName);
            FileReader fileReader=new FileReader(newFile);
            BufferedReader reader=new BufferedReader(fileReader);
            String line = null;
            while ((line = reader.readLine()) != null) {
             inputs.add(line);
            }
            
            List<String> mainData= new <String>ArrayList();
            
            for(String str:inputs){
                if(str.charAt(0)=='#'){
                    out.println(str);
                }
                else
                    mainData.add(str);
            }
            int totalPeaks=mainData.size();
            //Calculating Element Number to be used for the sorted Scores array
            int elementNumber;
            if(filterLevel==0)
                elementNumber=0;
            else
                elementNumber= (totalPeaks*filterLevel/100) -1;
            
            //Checking for errors
            if(elementNumber>totalPeaks || elementNumber<0){
                System.out.println("Invalid Percentile Value, enter a Value between 0 and 100");
                System.out.println("Total Peaks:");
                System.out.println(totalPeaks);
                System.out.println("Requested Peak Number:");
                System.out.println(elementNumber);
                return;
            }
            
            System.out.println("Returning all peaks with percentile Score above "+Integer.toString(filterLevel)+"%");
            System.out.println("Total Peaks:");
            System.out.println(totalPeaks);
            System.out.println("Requested Peak Number:");
            System.out.println(elementNumber);    
            
            List<Double> scores=new <Double>ArrayList();
            
            for(String data:mainData){
                StringTokenizer dataElements=new StringTokenizer(data);
                //Score is stored in the 8th Token parsing through first 7 Tokens
                for(int i=0;i<7;i++){
                    dataElements.nextToken();
                }
                double score=Double.parseDouble(dataElements.nextToken());
                scores.add(score);
            }
            
            List <Double>sortedScores=mergeSort(scores);
            
            Double cutOffScore=sortedScores.get(elementNumber);
            System.out.println("CutOff Score:");
            System.out.println(cutOffScore);
            System.out.println("Maximum Score:");
            System.out.println(sortedScores.get(sortedScores.size()-1));
            System.out.println("Minimum Score:");
            System.out.println(sortedScores.get(0));
            
            for(int i=0;i<scores.size();i++){
                if(scores.get(i)>cutOffScore){
                    out.println(mainData.get(i));
                }
            }
            
            
            out.close();
        
        }
                
        catch(Exception e)
        {
         e.printStackTrace();
        }    
       
    }
    
}
