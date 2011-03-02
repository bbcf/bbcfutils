package ch.epfl.bbcf.sqlite;

import java.io.File;
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


/**
 * A java class which regroup methods to access and build a sqlite 
 * database, formated for the BBCF (for qualitative and quantitative data),
 * separated in two categories {MANIPULATION} & {CONSTRUCTION}
 * @author Yohan Jarosz
 *
 */
public class SQLiteAccess {

	/**
	 * a sqlite connection to a database
	 */
	private Connection connection;
	/**
	 * key = chromosome name,so the table name
	 * value = a prepared statement
	 */
	private Map<String,PreparedStatement> insertStatement;

	/**
	 * current number of queries
	 */
	private int nbQueries;

	/**
	 * the number limit of queries to reach
	 * before committing 
	 */
	private int limitQueries;


	private SQLiteAccess(Connection connection){
		this.connection = connection;
		this.insertStatement = new HashMap<String, PreparedStatement>();
		this.limitQueries = 100000;
		this.nbQueries = 0;
	}
	private SQLiteAccess(Connection conn,int limitQueries){
		this.connection = conn;
		this.insertStatement = new HashMap<String, PreparedStatement>();
		this.limitQueries = limitQueries;
		this.nbQueries=0;
	}

	/**
	 * return a SQLiteAccess object from a sqlite database
	 * with the connection to the fdatabase built in
	 * DON'T FORGET TO CALL close() METHOD WHEN YOU FINISHED TO USE SQLiteAccess
	 * @param fullPath - the path oh the database
	 * @return - SQLiteAccess
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
	 * with the connection to the fdatabase built in
	 * DON'T FORGET TO CALL close() METHOD WHEN YOU FINISHED TO USE SQLiteAccess
	 * @param fullPath - the path oh the database
	 * @param limitQueriesSize - the limit of queries to reach 
	 * before doing an automatic commit
	 * useful when you build a database and have high number of insert to do -
	 * set to 100000 by default 
	 * @return - SQLiteAccess
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLiteAccess
	 */
	public static SQLiteAccess getConnectionWithDatabase(String fullPath,int limitQueriesSize)
	throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException{
		Class.forName("org.sqlite.JDBC").newInstance();
		Connection conn = DriverManager.getConnection("jdbc:sqlite:/"+fullPath);
		conn.setAutoCommit(false);
		return new SQLiteAccess(conn,limitQueriesSize);
	}

	/**
	 * the method to call when you finished to use
	 * SQLiteAccess
	 * @throws SQLException
	 */
	public void close() throws SQLException{
		this.connection.close();
	}
	/**
	 * commit
	 * @throws SQLException
	 */
	public void commit() throws SQLException{
		this.connection.commit();
	}
	/**
	 *  execute the query
	 * @param statement - the statement
	 * @param query - the query
	 * @return a Resultset
	 * @throws SQLException
	 */
	public static ResultSet getResultSet(Statement statement,String query) throws SQLException{
		return statement.executeQuery(query);
	}


	//MANIPULATION METHODS//

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



	//CONSTRUCTION METHODS//
	/**
	 * {CONSTRUCTION}
	 * create a new database
	 * @param type : must be qualitative or quantitative
	 */
	public void createNewDatabase(String type) throws SQLException{
		Statement stat;
		stat = this.connection.createStatement();
		stat.executeUpdate("drop table if exists attributes;");
		stat.executeUpdate("create table attributes (key text, value text);");
		PreparedStatement prep = this.connection.prepareStatement("insert into attributes values (?,?);");
		prep.setString(1, "datatype");
		prep.setString(2,type);
		prep.executeUpdate();
	}


	/**
	 * {CONSTRUCTION}
	 * create a new table for a chromosome for 
	 * QUALITATIVE data
	 * initialize a prepared statement which you can insert values in the table
	 * (start,end,score,name,strand,attributes)
	 * @param chromosome - the chromosome name
	 * @throws SQLException
	 */
	public void newChromosome_qual(String chromosome) throws SQLException {
		PreparedStatement stat = this.connection.prepareStatement("create table if not exists " +
				""+protect(chromosome)+" (start integer,end integer,score real,name text,strand integer,attributes text);");
		stat.execute();
		this.connection.commit();
		this.insertStatement.put(
				chromosome,this.connection.prepareStatement("insert into "+protect(chromosome)+" values (?,?,?,?,?,?);"));

	}

	/**
	 * {CONSTRUCTION}
	 * create a new table for a chromosome for 
	 * QUANTITATIVE data
	 * initialize prepared statement which you can insert values in the table (start,end,score)
	 * @param chromosome - the chromosome name
	 *  
	 * @throws SQLException
	 */
	public void newChromosome_quant(String chromosome) throws SQLException{
		PreparedStatement stat = this.connection.prepareStatement("create table if not exists "+protect(chromosome)+" (start integer,end integer,score real);");
		stat.execute();
		this.connection.commit();
		this.insertStatement.put(
				chromosome,this.connection.prepareStatement("insert into "+protect(chromosome)+" values (?,?,?);"));
	}


	/**
	 * {CONSTRUCTION}
	 * get the all the chromosomes names in a sqlite database before the table
	 * chromosome names is created
	 * @return - a list of chromosome names
	 * @throws SQLException 
	 */
	public List<String> getChromosomesNames() throws SQLException {
		List<String> chrNames = new ArrayList<String>();
		Statement stat = connection.createStatement();
		String query = "SELECT name FROM sqlite_master where type='table'and name!='attributes';";
		ResultSet rs = getResultSet(stat, query);
		while (rs.next()) {
			chrNames.add(rs.getString(1));
		}
		rs.close();
		return chrNames;
	}
	/**
	 * {CONSTRUCTION}
	 * look if the chromosome is already created
	 * @param chromosome - the chromosome name
	 * @return - true if the chromosome is created
	 */
	public boolean isCromosomeCreated(String chromosome){
		return insertStatement.get(chromosome)!=null;
	}
	/**
	 * {CONSTRUCTION}
	 * write values in the insert statement of the qualitative database
	 * YOU MUST CALL newChromosome_qual(String hromosome) before
	 * @param chr - the chromosome
	 * @param start
	 * @param end
	 * @param score
	 * @param name - name of the feature
	 * @param strand
	 * @throws SQLException
	 */
	public void writeValues_qual(String chr,int start,int end,float score,String name, int strand) throws SQLException{
		PreparedStatement prep = insertStatement.get(chr);	
		if(null!=prep){
			prep.setInt(1,start);
			prep.setInt(2,end);
			prep.setFloat(3,score);
			prep.setString(4, name);
			prep.setInt(5, strand);
			prep.execute();
			verifyNbQueries();
		} else {
			System.err.println(
					"the prepared statement for the chromosome "+chr+" is null, you must call newChromosome_qual("+chr+") before ");
		}
	}

	/**
	 * {CONSTRUCTION}
	 * write values in the insert statement of the qualitative database
	 * YOU MUST CALL newChromosome_quant(String chromosome) before
	 * @param chr
	 * @param start
	 * @param end
	 * @param score
	 * @throws SQLException
	 */
	public void writeValues_quant(String chr,int start,int end,float score) throws SQLException{
		PreparedStatement prep = insertStatement.get(chr);	
		if(null!=prep){
			prep.setInt(1,start);
			prep.setInt(2,end);
			prep.setFloat(3,score);
			prep.execute();
			verifyNbQueries();
		} else {
			System.err.println(
					"the prepared statement for the chromosome "+chr+" is null, you must call newChromosome_quant("+chr+") before ");
		}
	}

	/**
	 * {CONSTRUCTION}
	 * finalize the database with the creation of the
	 * table chrNames (name,length) and an index
	 * @param chr_length a Map with key = chromosome name and value = it's length
	 * @throws SQLException
	 */
	public void finalizeDatabase(Map<String,Integer> chr_length) throws SQLException{
		Statement stat = this.connection.createStatement();
		stat.executeUpdate("create table chrNames (name text, length integer);");
		
		PreparedStatement prep = this.connection.prepareStatement("insert into chrNames values (?,?);");
		int i=0;
		for(Map.Entry<String, Integer>entry : chr_length.entrySet()){
			i++;
			String chr = entry.getKey();
			Statement stat2 = this.connection.createStatement();
			stat2.executeUpdate("create index if not exist ind_name_"+i+" on "+chr+" (name);");
			int length = entry.getValue();
			prep.setString(1, chr);
			prep.setInt(2,length);
			prep.executeUpdate();
		}
		this.connection.commit();
	}


	//OTHERS//

	/**
	 * convenient method to protect a string for sqlite database
	 * e.g : if you have a table called 2micron, 
	 * sqlite will not create the table unless you protect it
	 * with double quotes
	 * @param str
	 * @return
	 */
	private static String protect(String str){
		return "\""+str+"\"";
	}

	/**
	 * Overrided method of Object
	 * in order to properly close the connection
	 * when garbage collector is called
	 */
	public void finalize(){
		try {
			close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		try {
			super.finalize();
		} catch (Throwable e) {
			e.printStackTrace();
		}
	}

	/**
	 * convenient method to commit if the number of
	 * queries during the connection
	 * exceed the limit
	 * @throws SQLException
	 */
	private void verifyNbQueries() throws SQLException{
		this.nbQueries++;
		if(nbQueries>limitQueries){
			this.commit();
			nbQueries = 0;
		}
	}


	//FOR THE DAEMON <tranform_to_sqlite>
	/**
	 * {CONSTRUCTION}
	 * convenient method which create a database for the use of 
	 * the transform_to_sqlite daemon
	 */
	public static boolean initJobsDatabase(File jobs) throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException{
		Class.forName("org.sqlite.JDBC").newInstance();
		Connection conn = DriverManager.getConnection("jdbc:sqlite:/"+jobs.getAbsolutePath());
		Statement stat = conn.createStatement();
		stat.execute("create table jobs " + 
				"(file text, " +
				"trackId text," + //the trackId in gdv database
				"tmpdir text," + //a tmp directory to delete after usage 
				"extension text," + //extension of the file
				"mail text," + // the mail to feedback (nomail if no feedback)
		"nrassemblyid text);");  // the assembly id of the genome in generep

		return true;
	}

}
