package ch.epfl.bbcf.bbcfutils.parser.exception;

public class ParsingException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3101202959880245742L;

	public ParsingException(final Exception e) {
		super(e);
	}

	public ParsingException(final Exception e, final String message, final int line) {
		super(e.getLocalizedMessage()+"  "+message+" at line "+line);
	}

	public ParsingException(final String message, final int line) {
		super(message+" at line "+line);
	}

}
