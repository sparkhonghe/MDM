/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Alone;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author HHU
 */
public class KmerCalculate_main {
    public static void main(String args[]) throws IOException {  
           // String inputFile = "D:"+File.separator+"AncestorSeq"+File.separator+"21"+File.separator+"21(33411239)38211238.fas";
            StringBuilder StrLine = new StringBuilder();   
            String outputFile = "D:"+File.separator+"Arabidopsis.txt";
            Output.toFile(outputFile, "The results come from KmerCalculate_main.java \n");
            HashMap<String, Long> KmersSetTotal = new HashMap<>();
            for(int j=5;j>0; j--){               
                //String inputFile = "reports"+File.separator+"refereces_all_chr"+File.separator+"fetchResults"+File.separator+j+File.separator+j+".fas"; //G:\MRAS\reports\refereces_all_chr\fetchResults\1
                String inputFile = "db1"+File.separator+"Arabidopsis"+File.separator+"Arabidopsis_chr"+j+".fasta";
                String StrTempBuff = "";
                int UnitSize = 5000;//以fa的行数为单位
                int n=0; //行计数
                int k=2;        
                HashMap<String , Long> KmersSet = new HashMap<>();
                File f = new File(inputFile);//读入.fas格式的参考序列
                LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
                String lineStr = lineNumberReader.readLine(); 
                //if(lineStr.startsWith(">")) lineStr = lineNumberReader.readLine(); 

                while(lineStr!=null){      
                    if(lineStr.startsWith(">")) {
                        lineStr = lineNumberReader.readLine();
                    }else{                      
                        StrTempBuff = StrTempBuff + lineStr;
                        n++;
                        if(n%UnitSize==0){
                            StrTempBuff = StrTempBuff.replaceAll("[\\s*\\t\\n\\r]", ""); //去掉.fas每行的回车符 
                            int seqLength = StrTempBuff.length();
                            if(seqLength>k){
                                KmerCalculate  newOne = new KmerCalculate(KmersSet);
                                KmersSet = newOne.getKmers(StrTempBuff,k);   
                                int newStartP = lineStr.length()-k+2;
                                StrTempBuff = lineStr.substring(newStartP, lineStr.length());//在两次reference sequence搜索时，为避免前一次处理和后一次处理中间衔接的地方有patternText存在                            
                            }
                        }
                    }
                    lineStr = lineNumberReader.readLine();  
                }

                if(StrTempBuff.length()!=0){
                        KmerCalculate  newOne = new KmerCalculate(KmersSet);
                        KmersSet = newOne.getKmers(StrTempBuff,k);                
                }
                System.out.println("CHR："+j);
                
                for(String  key:KmersSet.keySet()){
                    if(KmersSetTotal.get(key)!=null){
                        Long TempCount = Long.parseLong(KmersSetTotal.get(key).toString())+KmersSet.get(key);
                        KmersSetTotal.put(key, TempCount);
                    }else{
                            KmersSetTotal.put(key,KmersSet.get(key));
                    }
                    String TempStr = key+"\t"+KmersSet.get(key)+"\n";
                    StrLine.append(TempStr);
                    System.out.println(key+"\t"+KmersSet.get(key));
                }
                StrLine.append("\n\n");
                Output.append(outputFile, StrLine.toString());
                StrLine.setLength(0);                 
            }//染色体循环结束
            //System.out.println("通过Map.entrySet遍历key和value");
            for(Map.Entry<String,Long>entry:KmersSetTotal.entrySet()){
                String TempStr = entry.getKey()+"\t"+entry.getValue()+"\n";
                StrLine.append(TempStr);
                System.out.println(entry.getKey()+"\t"+entry.getValue());
            }
            Output.append(outputFile, StrLine.toString());
            StrLine.setLength(0);   
    }
    
}
