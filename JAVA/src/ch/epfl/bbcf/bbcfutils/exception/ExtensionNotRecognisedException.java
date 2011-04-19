package ch.epfl.bbcf.bbcfutils.exception;

import ch.epfl.bbcf.bbcfutils.conversion.sqlite.ConvertToSQLite.Extension;

public class ExtensionNotRecognisedException extends Exception{

	public ExtensionNotRecognisedException(Extension extension) {
		super(extension.name());
	}

}
