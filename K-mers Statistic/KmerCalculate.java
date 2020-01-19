/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Alone;

import java.util.HashMap;

/**
 *
 * @author HHU
 */
public class KmerCalculate {
    //String Sequence;
    //int k;
    HashMap<String , Long> KmersSet;
    public KmerCalculate(HashMap<String, Long> KmersSet) {
        this.KmersSet = KmersSet;
    }
    
    HashMap getKmers(String seq, int k){    
	int seqLength = seq.length();
	if(seqLength > k){
            for(int i = 0; i < seqLength - k + 1; i++){
                String TempStr = seq.substring(i, k + i).toUpperCase();
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
	}else {
		System.out.println(seq);
	}
        return KmersSet;
    }
}
