package ch.epfl.bbcf.bbcfutils.exception;

import ch.epfl.bbcf.bbcfutils.access.genrep.MethodNotFoundException;

public class ParsingException extends Exception{

	public ParsingException(NumberFormatException nfe, String string, int lineNb) {
		// TODO Auto-generated constructor stub
	}

	public ParsingException(String string, int lineNb) {
		// TODO Auto-generated constructor stub
	}

	public ParsingException(Exception e) {
		super(e.getMessage());
		
	}

	public ParsingException(String string) {
		// TODO Auto-generated constructor stub
	}

}
