package ch.epfl.bbcf.bbcfutils.access.gdv;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;


/**
 * 
 * @author Yohan Jarosz
 *
 */
public class PostToGDV {

	public static final String GDV_ADRESS = "http://svitsrv25.epfl.ch/gdv/post";
	private static final String GDV_POST = "gdv_post";

	/**
	 * Main entry of the program
	 * @param args
	 */
	public static void main(String[] args) {
		if(!parseArgs(args)){
			usage();
		}
	}




	/**
	 * parse the arguments provided by the user
	 * @param args - the args
	 * @param mail - the login
	 * @param pass - the password
	 * @return
	 */
	private static boolean parseArgs(String[] args) {
		String id = GDV_POST;
		if(args.length>0){
			if(args[0].equalsIgnoreCase("h")||args[0].equalsIgnoreCase("help")){
				System.out.println(README);	
				return false;
			} else if(args[0].equalsIgnoreCase("gendoc")){
				File file = new File("README");
				try {
					FileWriter writer = new FileWriter(file);
					writer.write(README);
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
				return true;
			}
		}
		RequestParameters params = buildRequestParameters(args);
		System.out.println("parsing args .... ");
		if(null!=params){
			if(checkParams(params.getMail(),params.getKey(),params.getCommand())){
				Command control = Command.getControl(params.getCommand(), params.getMail(), params.getKey());
				if(null!=control){
					return control.doRequest(id,params);
				}
			}
		}
		return false;
	}

	/**
	 * build an RequestParamters object with all args provided 
	 * by user
	 * @param args - the arguments
	 * @return RequestParamters
	 */
	private static RequestParameters buildRequestParameters(String[] args) {
		RequestParameters params = new RequestParameters();
		for(int i=0;i<args.length;i++){
			if(i+1<args.length){
				if(args[i].startsWith("--")){
					String cmd = args[i].substring(2, args[i].length());
					if(cmd.equalsIgnoreCase(RequestParameters.PROJECT_ID_PARAM)){
						params.setProjectId(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.URL_PARAM)){
						params.setUrl(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.COMMAND_PARAM)){
						params.setCommand(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.KEY_PARAM)){
						params.setKey(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.MAIL_PARAM)){
						params.setMail(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.DATATYPE_PARAM)){
						params.setDatatype(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.TYPE_PARAM)){
						params.setType(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.OBFUSCATED_PARAM)){
						params.setObfuscated(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.SEQUENCE_ID_PARAM)){
						params.setSequenceId(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.NAME_PARAM)){
						params.setName(args[i+1]);
					} else if(cmd.equalsIgnoreCase(RequestParameters.PUBLIC_PARAM)){
						params.setIsPublic(args[i+1]);
					} else {
						System.err.println("the parameter "+args[i]+" is not a parameter.");
					}
				} else {
					System.err.println(args[i]+" is not a parameter. Perhaps you didn't add --.");
				}
			} else {
				System.err.println("the parameters "+args[i]+" is not defined");
				return null;
			}
			i++;
		}
		return params;
	}







	/**
	 * check if a parameter required is given or not
	 * @param params
	 * @return
	 */
	public static boolean checkParams(String... params) {
		for(String p : params){
			if(null==p){
				System.err.println("\nmissing param(s)");
				return false;
			}
		}
		return true;
	}

	/**
	 * convenient method to write a table to the output
	 * @param values
	 * @return
	 */
	private static String writeList(String[] values) {
		String str = "";
		for(String v :values){
			str+=v+", ";
		}
		str=str.substring(0, str.length()-2);
		return str;
	}

	/**
	 * print help to the console output
	 */
	private static void usage() {
		String str = 
			"\n-----------" +
			"\n-----------" +
			"\n--SUMMARY--" +
			"\n-----------" +
			"\n-----------" +
			"\nthe args are given like;\n" +
			"\t --<arg> <value> (without <>) " +
			"e.g: --name Bob";
		System.out.println(str);
		System.out.println();
		System.out.println("--REQUESTED PARAMETERS : ");
		for(Map.Entry<String,String[]> entry : RequestParameters.getMapcommands().entrySet()){
			System.out.println(entry.getKey()+" : ");
			System.out.println("\t "+writeList(entry.getValue()));
		}
		System.out.println("-----------\n-----------");
	}

	private static final String README =
		"\n################" +
		"\n#### README ####" +
		"\n################" +
		"\n"+
		"\n#BUILD" +
		"\ndownload project, cd project directory and run the comment \"ant jar\"" +
		"\n"+
		"\n#INTRO" +
		"\nFirst you need a mail and an user key." +
		"\nLogin in GDV and go to 'preferences' page you will have your personal user key." +
		"\nNow there is two way : you can use the jar file provided or send by yourself" +
		"\na POST connection to \"http://svitsrv25.epfl.ch/gdv/post\". Note that if you send a POST" +
		"\nby yourself you will have to add the parameter \"id=gdv_post\" in your request." +
		"\n"+
		"\n#DEFINITION OF THE PARAMETERS" +
		"\n"+
		"\nREQUESTED :" +
		"\n\t- mail : login in tequila" +
		"\n\t- key  : user key" +
		"\n\t- command : one of [new_project, add_track,add_sqlite]" +
		"\n"+
		"\nREQUESTED BY COMMAND :" +
		"\n->'new_project' :" +
		"\n\t- seq_id : the nr_assembly id of the species in Genrep" +
		"\n\t- name : name of the project" +
		"\n\t- public (OPTIONNAL default=false): true to make the project public" +
		"\n\tRESPONSE : a JSON String : {\"project_id\":<anID>,\"public_url\":<anURL>}" +
		"\n"+
		"\n->'add_track' :" +
		"\n\t- url : url to fetch" +
		"\n\t- project_id : the project id in GDV" +
		"\n\t- name (OPTIONNAL default=name of the file): the name you want to give to the track" +
		"\n"+
		"\n->'add_sqlite' :" +
		"\n\t- url : url to fetch" +
		"\n\t- project_id : the project id in GDV" +
		"\n\t- datatype : one of [qualitative,quantitative] N.B :QUALITATIVE does not work yet" +
		"\n\t- name (OPTIONNAL default=name of the file): the name you want to give to the track" +
		"\n"+
		"\n#REQUESTS EXAMPLES" +
		"\n1) Want to create a public project in your profile for yeast?" +
		"\n\t$java -jar --mail <your mail> -key <your user key> --command new_project --seq_id 98 --name TestProject --public true" +
		"\n\t\t- or -"+
		"\n\tPOST at http://svitsrv25.epfl.ch/gdv_dev/post" +
		"\n\twith body = \"id=gdv_post&mail=<your mail>&key=<your user key>&command=new_projectseq_id=98&name=TestProject&public=true " +
		"\n\tIt will give you back an ID and an URL if you made it public. This is the project id in GDV and	the URL	you can	give to	" +
		"\n\tsomeone	else (e.g. {\"project_id\":5622,\"public_url\":\"http://svitsrv25.epfl.ch/gdv/browser?id=5622&ukey=ol9iio7k11&pkey=sg67v1ceqo83\"}" +
		"\n"+
		"\n2) Want to add a track to a project ?" +
		"\n\t$java -jar  --mail <your mail> -key <your user key> --command add_track --url <http://salt/BED/myfile.wig> --project_id <project id>" +
		"\n\tfor POST, see (1).  (don't forget to add id=gdv_post to the request's body)" +
		"\n"+
		"\n3) Want to add a track in sqlite format ?" +
		"\n\t$java -jar  --mail <your mail> -key <your user key> --command add_sqlite --url <http://salt/BED/myfile.db> --project_id <project id> --datatype quantitative" +
		"\n	(qualitative not working for the moment)";
}

