package ch.epfl.bbcf.bbcfutils.access.gdv;

import java.io.IOException;


import ch.epfl.bbcf.bbcfutils.access.InternetConnection;
import ch.epfl.bbcf.bbcfutils.access.gdv.RequestParameters.CMD;
import ch.epfl.bbcf.bbcfutils.access.gdv.RequestParameters.PARAM;






public class Command {

	public static boolean doRequest(String url,RequestParameters params) {
		/* prepare the body */
		String body = addParam(PARAM.id,params.getId(),true);
		body+=addParam(PARAM.command,params.getCommand().name());
		body+=addParam(PARAM.mail,params.getMail());
		body+=addParam(PARAM.key,params.getKey());


		CMD command = params.getCommand();
		boolean ok = true;
		switch(command){

		case new_project:/* launch a new project on GDV */
			System.out.println("NEW p");
			ok = checkParam(PARAM.name,params.getName());
			ok = checkParam(PARAM.seq_id,params.getSequenceId());
			body+=addParam(PARAM.name,params.getName());
			body+=addParam(PARAM.seq_id,params.getSequenceId());
			body+=addParam(PARAM.Public, params.getIsPublic());
			break;
			
			
		case new_track:/* add a track on a GDV project */
			ok = checkParam(PARAM.project_id,params.getProjectId());
			ok = checkParam(PARAM.url,params.getUrl());
			body+=addParam(PARAM.project_id,params.getProjectId());
			body+=addParam(PARAM.url,params.getUrl());
			body+=addParam(PARAM.name,params.getName());
			break;
			
			
		case assemblies:/* look at assemblies created in GDV */
			break;
			
			
		case status :/* get the status of a job */
			ok = checkParam(PARAM.job_id,params.getJob_id());
			body+=addParam(PARAM.job_id,params.getJob_id());
			break;
		default :
			System.err.println("command not recognized ");
		}
		if(ok){
			try {
				String result = InternetConnection.sendPOSTConnection(url, body, InternetConnection.MIME_TYPE_FORM_APPLICATION);
				System.out.println(result);
				return true;
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return false;
	}

	
	
	
	
	
	
	
	
	private static String addParam(PARAM key,String value,boolean first){
		String s="";
		if(!first){
			s+="&";
		}
		s+=key+"="+value;
		return s;
	}

	private static String addParam(PARAM key,String value){
		return addParam(key, value, false);
	}


	public static boolean checkParam(final PARAM p,final String value) {
		if(null==value){
			System.err.println("missing param : "+p.name());
			return false;
		}
		return true;
	}
	public static boolean checkCommand(final CMD value) {
		if(null==value){
			System.err.println("missing command");
			return false;
		}
		return true;
	}
}

