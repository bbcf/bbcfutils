package ch.epfl.bbcf.bbcfutils.conversion.json;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.json.JSONException;

import ch.epfl.bbcf.bbcfutils.conversion.json.ConvertToJSON.JSONFeat;



/**
 * class which sort JSONFeature, and arrange them 
 * like an NCList see http://bioinformatics.oxfordjournals.org/content/23/11/1386.full.pdf
 * @author Yohan Jarosz
 *
 */
public class NCList {
	public static void sort( List<JSONFeat> sortables){
		Collections.sort(sortables);
	}

	
	public static List<JSONFeat> arrange(List<JSONFeat> list) throws JSONException{
		List<JSONFeat> endList = new ArrayList<JSONFeat>();
		for(int i=0;i<list.size();i++){
			JSONFeat curFeat = list.get(i);
			i = includeNextFeat(list,i,curFeat,true);
			endList.add(curFeat);
		}
		return endList;
	}


	private static int includeNextFeat(List<JSONFeat> list, int i, JSONFeat curFeat,boolean first) throws JSONException {
		if(first){
			curFeat.initFeature();
		}
		if(i+1<list.size()){
			JSONFeat nextFeat = list.get(i+1);
			if(curFeat.getEnd()>nextFeat.getEnd()){
				i++;
				i = includeNextFeat(list,i,nextFeat,true);
				curFeat.nesting(nextFeat);
				i = includeNextFeat(list,i,curFeat,false);
			}
		}
		return i;
	}




}
