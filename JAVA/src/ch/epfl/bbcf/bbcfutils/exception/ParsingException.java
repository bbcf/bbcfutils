package ch.epfl.bbcf.bbcfutils.exception;

import ch.epfl.bbcf.bbcfutils.access.genrep.MethodNotFoundException;

public class ParsingException extends Exception{

	public ParsingException(NumberFormatException nfe, String string, int lineNb) {
		super(nfe.getMessage()+"   "+string+"  at line "+lineNb);
	}

	public ParsingException(String string, int lineNb) {
		super(string+"  at line "+lineNb);
	}

	public ParsingException(Exception e) {
		super(e.getMessage());
		
	}

	public ParsingException(String string) {
		super(string);
	}

}
