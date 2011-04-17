package ch.epfl.bbcf.bbcfutils.sqlite;

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

import ch.epfl.bbcf.bbcfutils.parser.feature.ExtendedQualitativeFeature;
import ch.epfl.bbcf.bbcfutils.parser.feature.QualitativeFeature;


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
	 * {MANIPULATION}
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
	 * {MANIPULATION}
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

	
	public ResultSet prepareQualitativeFeatures(String chr) throws SQLException {
		String query = "SELECT * FROM "+protect(chr)+" order by start asc";
		Statement prep = connection.createStatement();
		ResultSet r = prep.executeQuery(query);
		return r;
	}


	public QualitativeFeature getNextQualitativeFeature(ResultSet r) throws SQLException{
		QualitativeFeature feat = new QualitativeFeature();
		feat.setStart(r.getInt(1));
		feat.setEnd(r.getInt(2));
		feat.setScore(r.getFloat(3));
		feat.setName(r.getString(4));
		feat.setStrand(r.getInt(5));
		feat.setAttributes(r.getString(6));
		return feat;
	}
	public ExtendedQualitativeFeature getNextExtendedQualitativeFeature(ResultSet r) throws SQLException{
		ExtendedQualitativeFeature feat = new ExtendedQualitativeFeature();
		feat.setStart(r.getInt(1));
		feat.setEnd(r.getInt(2));
		feat.setScore(r.getFloat(3));
		feat.setName(r.getString(4));
		feat.setStrand(r.getInt(5));
		feat.setAttributes(r.getString(6));
		feat.setType(r.getString(7));
		feat.setIdentifier(r.getString(8));
		return feat;
	}

	

}
