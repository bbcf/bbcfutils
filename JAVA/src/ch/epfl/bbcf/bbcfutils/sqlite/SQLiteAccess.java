package ch.epfl.bbcf.bbcfutils.sqlite;

import java.io.BufferedReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ch.epfl.bbcf.bbcfutils.exception.ExtensionNotRecognisedException;
import ch.epfl.bbcf.bbcfutils.parsing.SQLiteExtension;
import ch.epfl.bbcf.bbcfutils.parsing.feature.BioSQLiteQualitative;
import ch.epfl.bbcf.bbcfutils.parsing.feature.BioSQLiteQualitativeExt;
import ch.epfl.bbcf.bbcfutils.parsing.feature.BioSQLiteQuantitative;


/**
 * A java class which regroup methods to access a sqlite 
 * database, formated for the BBCF (for qualitative and quantitative data),
 * @author Yohan Jarosz
 *
 */
public class SQLiteAccess extends SQLiteParent{



	protected SQLiteAccess(Connection connection) {
		super(connection);
	}


	protected SQLiteAccess(Connection conn, int limitQueries) {
		super(conn, limitQueries);
	}

	/**
	 * return a SQLiteAccess object from a sqlite database
	 * with the connection to the database built in
	 * DON'T FORGET TO CALL close() METHOD WHEN YOU FINISHED TO USE SQLiteAccess
	 * @param fullPath - the path oh the database
	 * @return - SQLiteConstruct
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLiteAccess
	 */
	public static SQLiteAccess getConnectionWithDatabase(String fullPath)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException{
		Class.forName("org.sqlite.JDBC").newInstance();
		Connection conn = DriverManager.getConnection("jdbc:sqlite:/"+fullPath);
		conn.setAutoCommit(false);
		return new SQLiteAccess(conn,100000);
	}
	/**
	 * return a SQLiteAccess object from a sqlite database
	 * with the connection to the database built in
	 * DON'T FORGET TO CALL close() METHOD WHEN YOU FINISHED TO USE SQLiteAccess
	 * @param fullPath - the path oh the database
	 * @param limitQueriesSize - the limit of queries to reach 
	 * before doing an automatic commit
	 * useful when you build a database and have high number of insert to do -
	 * set to 100000 by default 
	 * @return - SQLiteConstruct
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLiteAccess
	 */
	protected static SQLiteAccess getConnectionWithDatabase(String fullPath,int limitQueriesSize)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException{
		Class.forName("org.sqlite.JDBC").newInstance();
		Connection conn = DriverManager.getConnection("jdbc:sqlite:/"+fullPath);
		conn.setAutoCommit(false);
		return new SQLiteAccess(conn,limitQueriesSize);
	}


	/**
	 * get the maximum score for a chromosome
	 * @param chromosome
	 * @return a Float
	 * @throws SQLException
	 */
	public Float getMaxScoreForChr(String chromosome) throws SQLException{
		Statement stat = connection.createStatement();
		String query = "SELECT max(score) FROM "+protect(chromosome)+";";
		ResultSet rs = getResultSet(stat, query);
		Float f = null;
		while (rs.next()) {
			f= rs.getFloat(1);
		}
		rs.close();
		return f;
	}

	/**
	 * get the minimum score for a chromosome
	 * @param chromosome
	 * @return a Float
	 * @throws SQLException
	 */
	public float getMinScoreForChr(String chromosome) throws SQLException {
		Statement stat = connection.createStatement();
		String query = "SELECT min(score) FROM "+protect(chromosome)+";";
		ResultSet rs = getResultSet(stat, query);
		Float f = null;
		while (rs.next()) {
			f= rs.getFloat(1);
		}
		rs.close();
		return f;
	}


	/**
	 * {MANIPULATION}
	 * Take the chromosome & theirs length
	 * @return an Map with key : the chromosome name and value it's length
	 * @throws SQLException
	 */
	public Map<String, Integer> getChromosomesAndLength() throws SQLException {
		Map<String,Integer> result = new HashMap<String,Integer>();
		Statement stat = connection.createStatement();
		String query = "SELECT * FROM chrNames;";
		ResultSet rs = getResultSet(stat, query);
		while (rs.next()) {
			result.put(rs.getString("name"),
					rs.getInt("length"));
		}
		rs.close();
		return result;
	}

	/**
	 * get a value from the key you provide
	 * @param key
	 * @return a String
	 * @throws SQLException
	 */
	public String getAttribute(String key) throws SQLException {
		String result = null;
		String query = "select * from attributes where key = ? limit 1;";
		PreparedStatement prep = connection.prepareStatement(query);
		prep.setString(1, key);
		ResultSet rs = prep.executeQuery();
		while (rs.next()) {
			result = rs.getString("value");
		}
		rs.close();
		return result;
	}

	/**
	 * get all attributes of the file
	 * @return a map
	 * @throws SQLException
	 */
	public Map<String,String> getAttributes() throws SQLException {
		Map<String,String> atts = new HashMap<String, String>();
		String query = "select * from attributes;";
		PreparedStatement prep = connection.prepareStatement(query);
		ResultSet rs = prep.executeQuery();
		while (rs.next()) {
			atts.put(rs.getString(1), rs.getString(2));
		}
		rs.close();
		return atts;
	}

	public SQLiteExtension getDatatype() throws SQLException, ExtensionNotRecognisedException{
		String datatype = getAttribute("datatype");
		if(datatype.equalsIgnoreCase("qualitative")){
			return SQLiteExtension.QUALITATIVE;
		}
		if(datatype.equalsIgnoreCase("quantitative")){
			return SQLiteExtension.QUANTITATIVE;
		}if(datatype.equalsIgnoreCase("qualitative_extended")){
			return SQLiteExtension.QUALITATIVE_EXTENDED;
		}
		throw new ExtensionNotRecognisedException(datatype);
	}

	/**
	 * {MANIPULATION}
	 * Count the number of feature a chromosome have
	 * @param chromosome
	 * @return the count, or 0 if no feature found
	 * @throws SQLException 
	 */
	public int getFeatureCountForChromosome(String chromosome) throws SQLException {
		PreparedStatement prep = connection.prepareStatement("select count(*) from "+protect(chromosome)+";");
		ResultSet r = prep.executeQuery();
		if(r.next()){
			int result = r.getInt(1); 
			r.close();
			return result;
		}
		return 0;
	}

	/**
	 * {MANIPULATION}
	 * get the length of a chromosome
	 * in the database specified
	 * @param chr - the chromosome
	 * @return an int
	 * @throws SQLException 
	 */
	public int getLengthForChromosome(String chr) throws SQLException {
		PreparedStatement stat = connection.prepareStatement("SELECT length FROM chrNames where name = ? limit 1;");
		stat.setString(1, chr);
		ResultSet rs = stat.executeQuery();
		while (rs.next()) {
			int result = rs.getInt(1); 
			rs.close();
			return result;	
		}
		rs.close();
		return -1;
	}

	/**
	 * {MANIPULATION}
	 * get the maximum end of the chromosome's features
	 * in the database specified
	 * @param chromosome
	 * @return an int
	 * @throws SQLException 
	 */
	public int getMaxEndForChromosome(String chromosome) throws SQLException {
		int f = 0;
		Statement stat = connection.createStatement();
		String query = "SELECT max(end) FROM "+protect(chromosome)+";";
		ResultSet rs = getResultSet(stat, query);
		while (rs.next()) {
			f= rs.getInt(1);
		}
		rs.close();
		return f;
	}

	/**
	 * {MANIPULATION}
	 * get all start & end position for the chromosome
	 * @param chr - the chromosome
	 * @return a sql resltset that you have to close after use
	 * @throws SQLException
	 */
	public ResultSet getStartEndForChromosome(String chr) throws SQLException {
		Statement stat = connection.createStatement();
		String query = "SELECT start,end FROM \""+chr+"\"; ";
		ResultSet r = getResultSet(stat, query);
		return r;
	}

	/**
	 * try to find coordinates (start,end) of a gene by it's name
	 * - EXACT MATCH -
	 * @param chr - the chromosome
	 * @param name - the gene name
	 * @return a list of coordinates (start,end,start,end,start,end,....)
	 * @throws SQLException
	 */
	public List<Integer> searchForGeneNameOnChromosome(String chr,String name) throws SQLException {
		List<Integer> result = new ArrayList<Integer>();
		String query = "SELECT start,end FROM "+protect(chr)+" where name = ?; ";
		PreparedStatement prep = connection.prepareStatement(query);
		prep.setString(1, name);
		ResultSet r = getResultSet(prep, query);
		while(r.next()){
			result.add(r.getInt(1));
			result.add(r.getInt(2));
		}
		r.close();
		return result;
	}



	public SQLiteExtension getDatabaseDatatype() throws SQLException{
		String q1 = "SELECT * FROM attributes where key = 'datatype' limit 1; ";
		Statement stat = connection.createStatement();
		ResultSet r = stat.executeQuery(q1);
		if(r.next()){
			String dt = r.getString("value");
			return SQLiteExtension.valueOf(dt.toUpperCase());
		}
		return null;
	}



	public ResultSet prepareFeatures(String chr) throws SQLException {
		String query = "SELECT * FROM "+protect(chr)+" order by start asc";
		Statement prep = connection.createStatement();
		ResultSet r = prep.executeQuery(query);
		return r;
	}


	public BioSQLiteQuantitative getNextQuantitativeFeature(ResultSet r, String chr) throws SQLException{
		BioSQLiteQuantitative feat = new BioSQLiteQuantitative();
		feat.setChromosome(chr);
		feat.setStart(r.getInt(1));
		feat.setEnd(r.getInt(2));
		feat.setScore(r.getFloat(3));
		return feat;
	}

	public BioSQLiteQualitative getNextQualitativeFeature(ResultSet r, String chr) throws SQLException {
		BioSQLiteQualitative feat = new BioSQLiteQualitative();
		feat.setChromosome(chr);
		feat.setStart(r.getInt(1));
		feat.setEnd(r.getInt(2));
		try {
			feat.setScore(r.getFloat(3));
		} catch (SQLException e) {}
		try {
			feat.setName(r.getString(4));
		} catch (SQLException e) {}
		try {
			feat.setStrand(r.getInt(5));
		} catch (SQLException e) {}
		try {
			feat.setAttributes(r.getString(6));
		} catch (SQLException e) {}
		return feat;
	}
	public BioSQLiteQualitativeExt getNextExtendedQualitativeFeature(ResultSet r,String chr) throws SQLException{
		BioSQLiteQualitativeExt feat = new BioSQLiteQualitativeExt();
		feat.setChromosome(chr);
		feat.setStart(r.getInt(1));
		feat.setEnd(r.getInt(2));
		try {
			feat.setScore(r.getFloat(3));
		} catch (SQLException e) {}
		try {
			feat.setName(r.getString(4));
		} catch (SQLException e) {}
		try {
			feat.setStrand(r.getInt(5));
		} catch (SQLException e) {}
		try {
			feat.setAttributes(r.getString(6));
		} catch (SQLException e) {}
		try {
			feat.setType(getTypeForExtendedQualitativeFeature(
					r.getInt(7)));
		} catch (SQLException e) {}
		try {
			feat.setIdentifier(r.getString(8));
		} catch (SQLException e) {}
		return feat;
	}

	public BioSQLiteQualitativeExt getExtendedQualitativeFeature(String chr) throws SQLException{
		BioSQLiteQualitativeExt feat = new BioSQLiteQualitativeExt();
		String query = "SELECT * FROM "+protect(chr)+" limit 1";
		Statement prep = connection.createStatement();
		ResultSet r = prep.executeQuery(query);
		while(r.next()){
			feat.setStart(r.getInt(1));
			feat.setEnd(r.getInt(2));
			feat.setScore(r.getFloat(3));
			feat.setName(r.getString(4));
			feat.setStrand(r.getInt(5));
			feat.setAttributes(r.getString(6));
			feat.setType(
					getTypeForExtendedQualitativeFeature(r.getInt(7)));
			feat.setIdentifier(r.getString(8));
		}
		return feat;
	}

	public String getTypeForExtendedQualitativeFeature(int type) throws SQLException{
		String query = "select type from types where identifier = ? limit 1";
		PreparedStatement prep = connection.prepareStatement(query);
		prep.setInt(1, type);
		ResultSet r = prep.executeQuery();
		String t="default";
		if(r.next()){
			t=r.getString(1);
		}
		r.close();
		return t;
	}
	public List<String> getTypes() throws SQLException{
		String query = "select type from types";
		PreparedStatement prep = connection.prepareStatement(query);
		ResultSet r = prep.executeQuery();
		List<String> types = new ArrayList<String>();
		while(r.next()){
			types.add(r.getString(1));
		}
		r.close();
		return types;
	}




}
