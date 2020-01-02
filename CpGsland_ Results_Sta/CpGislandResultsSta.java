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

/**
 *
 * @author HHU
 */
public class CpGislandResultsSta {
    String inputFolder;
    
    public CpGislandResultsSta(String inputFolder) {
        this.inputFolder = inputFolder;
    }
    
    public Long[] StaMotif() throws IOException {//
        Long[] PatternNum = new Long[16];
        for(int i=0; i<PatternNum.length; i++) PatternNum[i]=0L;
        
        File inputFiles = new File(inputFolder);
        String[] filelist = inputFiles.list();
        for(int i=0; i<filelist.length; i++){
            String inputFileTemp = inputFolder+File.separator+filelist[i];
            File f = new File(inputFileTemp);
            LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(f));
            String lineStr = lineNumberReader.readLine();
            while (lineStr != null) { 
                String[] Token = lineStr.split("\\s+");
                Long NumTemp = 0L;
                for(int j=9; j<Token.length; j++){
                    PatternNum[j-9] = PatternNum[j-9] + Long.parseLong(Token[j]);
                    NumTemp = NumTemp + Long.parseLong(Token[j]);
                }
                if(NumTemp!=Long.parseLong(Token[5]))  System.out.println("Becare for!"+"\t"+NumTemp+"\t"+Long.parseLong(Token[5]));
                lineStr = lineNumberReader.readLine();
            }
        }
        return PatternNum;
    }
}
