/*
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
  * To change this license header, choose License Headers in Project Properties.
*/

package Alone;

/**
 *
 * @author MZ
 */


import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Motif_Finder {
    String inputFile;
    String outputFile;
    String chr;
    String patternText;
    String PatternList;

    public Motif_Finder(String inputFile, String outputFile, String chr) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
        this.chr = chr;        
    }
    public Motif_Finder(String inputFile, String outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
    }
/*
    public Motif_Finder(String inputFile, String patternText) { //for run2()
        this.inputFile = inputFile;
        this.patternText = patternText;
    }
  //  int Size = 1000000;
    */
    public Motif_Finder(String inputFile, String outputFile,String patternText,String chr) { //for run()
        this.inputFile = inputFile;
        this.outputFile = outputFile;
        this.patternText = patternText;
        this.chr = chr;
    }
    public int[][] CpGIsland_group(String CpGIslandFile) throws IOException {//建议输入文件为每条染色体的参考序列文件（没划分成子文件的），本程序将按给字的patternText进行搜索，找出参考序列中所包含的所有patternText的个数 
        List<String> CpGIslands = new ArrayList<String>();
        
        File f = new File(CpGIslandFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
            //System.out.println(lineStr);
            String[] Token = lineStr.split("\\s+");
            String Rec = Token[1]+"\t"+Token[2];
            CpGIslands.add(Rec);       
            lineStr = lineNumberReader.readLine();            
        }//****************************
        int[][] CpGArray =new int[CpGIslands.size()][3];
        for(int i=0; i<CpGIslands.size(); i++){
            String[] TempToken = CpGIslands.get(i).split("\\s+");
            CpGArray[i][0] = Integer.parseInt(TempToken[0]);
            CpGArray[i][1] = Integer.parseInt(TempToken[1]);   
            CpGArray[i][2] = 0; 
        }
        return CpGArray;
    } 
    public boolean CpGIslandJudge(int Position, int[][]CpGArray) throws IOException {//建议输入文件为每条染色体的参考序列文件（没划分成子文件的），本程序将按给字的patternText进行搜索，找出参考序列中所包含的所有patternText的个数 
        boolean rel = false;
        int Size = CpGArray.length;
        for(int i=0; i<Size; i++){
            if(Position>CpGArray[i][0]&&Position<CpGArray[i][1]) {
                rel = true;
                CpGArray[i][2]++;
                //System.out.println( CpGArray[i][0]+"\t"+ CpGArray[i][1]+"\t"+ CpGArray[i][2]);                
                return rel;
            }
        }
        return rel;
    }     
    public void FinderMotif_group_AvoidCpGIsland(String[] patternTextArray,int[][] CpGArray) throws IOException {//建议输入文件为每条染色体的参考序列文件（没划分成子文件的），本程序将按给字的patternText进行搜索，找出参考序列中所包含的所有patternText的个数 
        int[] NonIslandCountArray = new int[patternTextArray.length];
        int[] IslandCountArray = new int[patternTextArray.length];
        for(int i=0; i<patternTextArray.length; i++){
            NonIslandCountArray[i]=0;
            IslandCountArray[i]=0;
        }
        File f = new File(inputFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();
        int DealSize = 30000;
        StringBuilder StrLine = new StringBuilder();   
        int n=0;//行数计数器
        int lineNumb =0;
        int CGNumber =0;//CG含量计数器
        int IslandCGNumber =0;
        String TempStr = "";
        int StartPosition=0;
        int EndPosition = 0;
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
            if(lineStr.startsWith(">")){

            }else{
                //System.out.println("We are here!");
                TempStr = TempStr + lineStr;
                EndPosition = EndPosition + lineStr.length();
                n++; 
                if(n%DealSize==0){
                    //System.out.println("n="+n);
                    //System.out.println("Startposition="+StartPosition);
                    //System.out.println("EndPosition="+EndPosition);
                    TempStr = TempStr.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格
                    int PatternLength =0;
                    ///////////////////////////////////////////////////////////////
                    for(int i=0; i<NonIslandCountArray.length; i++){
                        CGNumber =0;//CG含量计数器
                        String patternTextTemp = patternTextArray[i];
                        Pattern pattern = Pattern.compile(patternTextTemp, Pattern.CASE_INSENSITIVE);
                        Matcher matcher = pattern.matcher(TempStr);

                        // using Matcher find(), group(), start() and end() methods
                        while (matcher.find()) {
                            int TempPosition = StartPosition + matcher.start();
                            if(CpGIslandJudge(TempPosition,CpGArray)){                                
                                IslandCGNumber++;
                            }else{
                                CGNumber++;
                            }
                        }
                        NonIslandCountArray[i] = NonIslandCountArray[i]+CGNumber;
                        IslandCountArray[i] = IslandCountArray[i]+IslandCGNumber;
                    }
                    ///////////////////////////////////////////////////////////
                   // System.out.println(patternTextArray[0]);
                    int lengthTemp = patternTextArray[0].length();
                    int newStartP = lineStr.length()-lengthTemp+2;
                    TempStr = lineStr.substring(newStartP, lineStr.length());//在两次reference sequence搜索时，为避免前一次处理和后一次处理中间衔接的地方有patternText存在
                    StartPosition = EndPosition -lengthTemp+2;
                   // System.out.println("StartPosition="+StartPosition);
                }
            }            
            lineStr = lineNumberReader.readLine();            
        }//****************************
        //不足处理数目的部分
        TempStr = TempStr.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格        
        ///////////////////////////////////////////////////////////////
        for(int i=0; i<NonIslandCountArray.length; i++){
            CGNumber =0;//CG含量计数器
            String patternTextTemp = patternTextArray[i];
            Pattern pattern = Pattern.compile(patternTextTemp, Pattern.CASE_INSENSITIVE);
            Matcher matcher = pattern.matcher(TempStr);

            // using Matcher find(), group(), start() and end() methods
            while (matcher.find()) {
                int TempPosition = StartPosition + matcher.start();
                if(!CpGIslandJudge(TempPosition,CpGArray)){
                    CGNumber++;
                }else{
                    IslandCGNumber++;
                }
            }
            NonIslandCountArray[i] = NonIslandCountArray[i]+CGNumber;
        }
        ///////////////////////////////////////////////////////////
        
        for(int i=0; i<NonIslandCountArray.length; i++){
            String newOne = patternTextArray[i]+"\t"+NonIslandCountArray[i]+"\t"+IslandCountArray[i]+"\n";
            StrLine.append(newOne);
        }
        Output.append(outputFile, StrLine.toString());
        StrLine.setLength(0);
        System.out.println("Total number of "+ patternText+ " is\t"+CGNumber + "\t in nonCpGIsland and\t"+IslandCGNumber+"\t in CpGIsland of chr "+ chr);
        
    } 
  
    public void FinderMotif_withoutIsland(String patternText,int[][] CpGArray) throws IOException {//建议输入文件为每条染色体的参考序列文件（没划分成子文件的），本程序将按给字的patternText进行搜索，找出参考序列中所包含的所有patternText的个数  
        File f = new File(inputFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();
        int DealSize = 1000;
        StringBuilder StrLine = new StringBuilder();   
        int n=0;//行数计数器
        int lineNumb =0;
        int CGNumber =0;//CG含量计数器
        String TempStr = "";
        int StartPosition=0;
        int FirstP=0;
        int EndP=0;
        int IslandCGNumber =0;
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
            if(lineStr.startsWith(">")){
                String[] Token = lineStr.split(":");
                StartPosition = Integer.parseInt(Token[3]);//取得该.fas的起始position，前提是文件的头行对应的地址一定要是正确的才行！！！！！！！！！！！！！！
                chr = Token[2];
            }else{
                //System.out.println("We are here!");
                TempStr = TempStr + lineStr;
                n++; 
                EndP = EndP + lineStr.length();
                if(n%DealSize==0){
                    TempStr = TempStr.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格
                    Pattern pattern = Pattern.compile(patternText, Pattern.CASE_INSENSITIVE);
                    Matcher matcher = pattern.matcher(TempStr);

                    // using Matcher find(), group(), start() and end() methods
                    while (matcher.find()) {
                        int TempPosition = StartPosition + matcher.start();
                        if(CpGIslandJudge(TempPosition,CpGArray)){                                
                            IslandCGNumber++;
                        }else{
                            CGNumber++;
                        }
                        //System.out.println("Found the text \"" + matcher.group() + "\" starting at " + matcher.start() + " index and ending at index " + matcher.end());
                        String subStrTemp = TempStr.substring(matcher.start()-1, matcher.end()+1);
                        String StrTemp = chr+"\t"+Integer.toString(matcher.start()+ StartPosition)+"\t"+Integer.toString(matcher.end()+ StartPosition)+"\t"+subStrTemp+"\n";
                        StrLine.append(StrTemp);
                        lineNumb++;
                        if(lineNumb%DealSize==0){
                            Output.append(outputFile, StrLine.toString());
                            StrLine.setLength(0);
                        }
                    }
                    int newStartP = lineStr.length()-patternText.length()-1;
                    TempStr = lineStr.substring(newStartP, lineStr.length());//在两次reference sequence搜索时，为避免前一次处理和后一次处理中间衔接的地方有patternText存在
                    FirstP = EndP -patternText.length()-1;
                    System.out.println(TempStr+"\t"+FirstP);
                }
            }            
            lineStr = lineNumberReader.readLine();            
        }//****************************
        //不足处理数目的部分
        TempStr = TempStr.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格
        Pattern pattern = Pattern.compile(patternText, Pattern.CASE_INSENSITIVE);
        Matcher matcher = pattern.matcher(TempStr);

        // using Matcher find(), group(), start() and end() methods
        while (matcher.find()) {
            int TempPosition = StartPosition + matcher.start();
            if(CpGIslandJudge(TempPosition,CpGArray)){                                
                IslandCGNumber++;
            }else{
                CGNumber++;
            }
            String subStrTemp = TempStr.substring(matcher.start()-1, matcher.end()+1);
            //System.out.println("Found the text \"" + matcher.group() + "\" starting at " + matcher.start() + " index and ending at index " + matcher.end());
            String StrTemp = chr+"\t"+Integer.toString(matcher.start()+ StartPosition)+"\t"+Integer.toString(matcher.end()+ StartPosition)+"\t"+subStrTemp+"\n";
            StrLine.append(StrTemp);
            lineNumb++;
            Output.append(outputFile, StrLine.toString());
            StrLine.setLength(0);
        }
        System.out.println("Total number of "+ patternText+ " is"+"\t"+CGNumber+"\tin non CpGIsland\t" +IslandCGNumber+ "\tin CpGIsland in chr "+ chr);
    }    
  
    public void FinderMotif() throws IOException {//建议输入文件为每条染色体的参考序列文件（没划分成子文件的），本程序将按给字的patternText进行搜索，找出参考序列中所包含的所有patternText的个数  
        File f = new File(inputFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();
        int DealSize = 10000;
        StringBuilder StrLine = new StringBuilder();   
        int n=0;//行数计数器
        int lineNumb =0;
        int CGNumber =0;//CG含量计数器
        String TempStr = "";
        int StartPosition=0;
        int EndPosition =0;
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
            if(lineStr.startsWith(">")){
                String[] Token = lineStr.split(":");
                StartPosition = Integer.parseInt(Token[4]);//取得该.fas的起始position，前提是文件的头行对应的地址一定要是正确的才行！！！！！！！！！！！！！！
                System.out.println("StartPosition = "+ StartPosition);
                EndPosition = StartPosition-1;
                chr = Token[2];
            }else{
                TempStr = TempStr + lineStr;
                n++; 
                EndPosition = EndPosition + lineStr.length();
                if(n%DealSize==0){
                    //System.out.println("StartPosition = "+ StartPosition);
                    TempStr = TempStr.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格
                    Pattern pattern = Pattern.compile(patternText, Pattern.CASE_INSENSITIVE);
                    Matcher matcher = pattern.matcher(TempStr);

                    // using Matcher find(), group(), start() and end() methods
                    while (matcher.find()) {
                        CGNumber++;
                        //System.out.println("Found the text \"" + matcher.group() + "\" starting at " + matcher.start() + " index and ending at index " + matcher.end());
                        String StrTemp = chr+"\t"+Integer.toString(matcher.start()+ StartPosition)+"\t"+Integer.toString(matcher.end()+ StartPosition)+"\n";
                        StrLine.append(StrTemp);
                        lineNumb++;
                        if(lineNumb%DealSize==0){
                            Output.append(outputFile, StrLine.toString());
                            StrLine.setLength(0);
                        }
                    }
                    int newStartP = lineStr.length()-patternText.length()+2;
                    TempStr = lineStr.substring(newStartP, lineStr.length());//在两次reference sequence搜索时，为避免前一次处理和后一次处理中间衔接的地方有patternText存在
                   // System.out.println("TempStr="+TempStr);
                   // System.out.println("StartPosition="+StartPosition);
                  //  System.out.println("EndPosition="+EndPosition);
                    StartPosition = EndPosition - TempStr.length();
                  //  System.out.println("NewStartPosition="+StartPosition);
                    EndPosition = StartPosition; 
                }
            }            
            lineStr = lineNumberReader.readLine();            
        }//****************************
        //不足处理数目的部分
        TempStr = TempStr.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格
        Pattern pattern = Pattern.compile(patternText, Pattern.CASE_INSENSITIVE);
        Matcher matcher = pattern.matcher(TempStr);

        // using Matcher find(), group(), start() and end() methods
        while (matcher.find()) {
            CGNumber++;
            //System.out.println("Found the text \"" + matcher.group() + "\" starting at " + matcher.start() + " index and ending at index " + matcher.end());
            String StrTemp = chr+"\t"+Integer.toString(matcher.start()+ StartPosition)+"\t"+Integer.toString(matcher.end()+ StartPosition)+"\n";
            StrLine.append(StrTemp);
            lineNumb++;
            Output.append(outputFile, StrLine.toString());
            StrLine.setLength(0);
        }
        System.out.println("Total number of "+ patternText+ " is "+CGNumber + "in chr "+ chr);
    }
    
    public Long[] FinderMotif_group(String[] patternTextArray) throws IOException {//建议输入文件为每条染色体的参考序列文件（没划分成子文件的），本程序将按给字的patternText进行搜索，找出参考序列中所包含的所有patternText的个数 
        HashMap<String, Long> KmersSet = new HashMap();
        Long[] patternNum = new Long[patternTextArray.length];
        for(int i=0; i<16;i++)  patternNum[i]=0L;
        
        File f = new File(inputFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();
        int DealSize = 10000;
          
        int n=0;//行数计数器
        String seq = "";
        int StartPosition=0;
	int EndPosition =0;
        int PatternLength = patternTextArray[0].length();
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
                if(lineStr.startsWith(">")){
                }else{
                    //System.out.println("We are here!");
                    seq = seq + lineStr;
                    n++; 
                    if(n%DealSize==0){
                            seq = seq.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格
                            ///////////////////////////////////////////////////////////////
                            int seqLength = seq.length();
                            if(seqLength > PatternLength){
                                for(int j = 0; j < seqLength - PatternLength + 1; j++){
                                    String TempStr = seq.substring(j, PatternLength + j).toUpperCase();
                                    if(KmersSet.get(TempStr)!=null){
                                        Long TempCount = Long.parseLong(KmersSet.get(TempStr).toString())+1;
                                        KmersSet.put(TempStr, TempCount);
                                    }else{
                                            if(!TempStr.contains("-")&&!TempStr.contains(".")&&!TempStr.contains("N")&&!TempStr.contains("n")){
                                                KmersSet.put(TempStr,1L);
                                            }else{
                                                //System.out.println(TempStr);
                                            }
                                    }
                                }                        

                            }
                            ///////////////////////////////////////////////////////////
                            int newStartP = lineStr.length()-PatternLength+2;
                            seq = lineStr.substring(newStartP, lineStr.length());//在两次reference sequence搜索时，为避免前一次处理和后一次处理中间衔接的地方有patternText存在
                           // System.out.println("TempStr="+TempStr);
                           // System.out.println("StartPosition="+StartPosition);
                          //  System.out.println("EndPosition="+EndPosition);                    
                    }            
           
                }//while end ****************************
                lineStr = lineNumberReader.readLine();  
        }
            //不足处理数目的部分
            //不足处理数目的部分
            seq = seq.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格     
            ///////////////////////////////////////////////////////////////
            int seqLength = seq.length();
            if(seqLength > PatternLength){
                for(int j = 0; j < seqLength - PatternLength + 1; j++){
                    String TempStr = seq.substring(j, PatternLength + j).toUpperCase();
                    if(KmersSet.get(TempStr)!=null){
                        Long TempCount = Long.parseLong(KmersSet.get(TempStr).toString())+1;
                        KmersSet.put(TempStr, TempCount);
                    }else{
                        if(!TempStr.contains("-")&&!TempStr.contains(".")&&!TempStr.contains("N")&&!TempStr.contains("n")){
                            KmersSet.put(TempStr,1L);
                        }else{
                            //System.out.println(TempStr);
                        }
                    }
                }                        
            }
            
            Long[] StrLine = new Long[patternTextArray.length];              
            for(int i=0; i<patternTextArray.length; i++){
                if(KmersSet.get(patternTextArray[i])!=null){
                    Long TempCount = Long.parseLong(KmersSet.get(patternTextArray[i]).toString());
                    //System.out.println(patternTextArray[i]+"\t"+TempCount);
                    String TempStr = patternTextArray[i]+"\t"+TempCount;
                    StrLine[i]=TempCount;
                }else{
                    //System.out.println(patternTextArray[i]+"\t"+"0");
                    String TempStr = patternTextArray[i]+"\t"+"0";
                    StrLine[i]=0L;
                }
            }
            return StrLine;
    }
 
        
    public void CG_ContainsContext() throws IOException {//本方法主要的目的是从Context_Get.java的输出结果文件中，搜索带有CpG-> TpG和CpG->CpA的所有Context,当然，也可以是ApT之类的二核苷
        String Non_CG = outputFile.replace(patternText, "Non_"+patternText);
        Output.toFile(Non_CG, "");
        File f = new File(inputFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();

        StringBuilder StrLine = new StringBuilder();  
        StringBuilder StrLine2 = new StringBuilder();
        int writeSize = 10000;
        int PatternTextNumber =0;//CG含量计数器
        int non_PatternTextNumber =0;
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
            String[] Token = lineStr.split("\\s+");   
            System.out.println("lineStr="+lineStr);
            int HalfLength = Token[0].length()/2;
            String newStr1 = Token[0].substring(HalfLength, HalfLength+2);
            String newStr2 = Token[0].substring(HalfLength-1, HalfLength+1);
            System.out.println("newStr1="+newStr1);
            System.out.println("newStr2="+newStr2);
            if(newStr1.equalsIgnoreCase(patternText)||newStr2.equalsIgnoreCase(patternText)) {
                StrLine.append(lineStr);
                StrLine.append("\n");
                PatternTextNumber++;
                if(PatternTextNumber%writeSize==0){
                    Output.append(outputFile, StrLine.toString());
                    StrLine.setLength(0);
                } 
            }else{
                StrLine2.append(lineStr);
                StrLine2.append("\n");  
                non_PatternTextNumber++;
                if(non_PatternTextNumber%writeSize==0){
                    Output.append(Non_CG, StrLine2.toString());
                    StrLine2.setLength(0);
                }                 
            }

            lineStr = lineNumberReader.readLine();            
        }//****************************
            Output.append(outputFile, StrLine.toString());
            StrLine.setLength(0); 
            
            Output.append(Non_CG, StrLine2.toString());
            StrLine2.setLength(0);
    }  
    
    public void run() throws IOException {//只处理二SNP，有两个或两个以上的突变等位基因的就会被跳过  
        File f = new File(inputFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();

        StringBuilder StrLine = new StringBuilder();   
        int n=0;//行数计数器
        int CGNumber =0;//CG含量计数器
        String TempStr = "";
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
            if(lineStr.startsWith(">")){
                //如果发现是.fas的起始行，就什么都不处理，略过
            }else{
                //System.out.println("We are here!");
                TempStr = TempStr + lineStr;
                n++; 
            }
            lineStr = lineNumberReader.readLine();            
        }//****************************

        TempStr = TempStr.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉行中的空格
        Pattern pattern = Pattern.compile(patternText, Pattern.CASE_INSENSITIVE);
        Matcher matcher = pattern.matcher(TempStr);

        // using Matcher find(), group(), start() and end() methods
        while (matcher.find()) {
            CGNumber++;
            System.out.println("Found the text \"" + matcher.group() + "\" starting at " + matcher.start() + " index and ending at index " + matcher.end());
        }  
        System.out.println("Total number of"+ patternText+ "is"+CGNumber);
    }
    
    public void run2() throws IOException {//本方法主要的目的是从Context_Get.java的输出结果文件中，搜索CpG-> TpG的突变，并计算这种突变对N->T的贡需比例 
        File f = new File(inputFile);
        LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
        String lineStr = lineNumberReader.readLine();

        StringBuilder StrLine = new StringBuilder();   
        int n=0;//行数计数器
        int cGTNumb = 0;
        int TNumb = 0;
        int CGNumber =0;//CG含量计数器
        String TempStr = "";
       // boolean C_onTail =  false;
        while (lineStr != null) { //**********************
            if(lineStr.startsWith(">")){
                //如果发现是.fas的起始行，就什么都不处理，略过
            }else{
                String[] Token = lineStr.split("\\s+"); 
                if(Token[1].equals("T")) TNumb = TNumb + Integer.parseInt(Token[3]);
                //System.out.println("We are here!");
                Pattern pattern = Pattern.compile(patternText);
                Matcher matcher = pattern.matcher(Token[0]);
                while (matcher.find()) {
                    CGNumber = CGNumber + Integer.parseInt(Token[3]);
                    if(Token[1].equals("T")){
                        cGTNumb = cGTNumb + Integer.parseInt(Token[3]);
                    }
                }    
                n++; 
            }
            lineStr = lineNumberReader.readLine();            
        }//****************************
        System.out.println("Total number of "+ patternText+ " is"+CGNumber);
        System.out.println("Total number of "+ patternText+ "change to TpG is"+cGTNumb);
        System.out.println("Total number of variations with Alter Allele T is "+TNumb);
    }
}

