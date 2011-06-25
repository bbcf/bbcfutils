package ch.epfl.bbcf.bbcfutils.exception;

import ch.epfl.bbcf.bbcfutils.parsing.Extension;


public class ExtensionNotRecognisedException extends Exception{

	public ExtensionNotRecognisedException(Extension extension) {
		super(extension.name());
	}

	public ExtensionNotRecognisedException(String datatype) {
		super(datatype);
	}

}
