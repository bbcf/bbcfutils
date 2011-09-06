package ch.epfl.bbcf.bbcfutils.sqlite;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import ch.epfl.bbcf.bbcfutils.parsing.SQLiteExtension;

/**
 * A java class which regroup methods to build a sqlite 
 * database, formated for the BBCF (for qualitative and quantitative data),
 * @author Yohan Jarosz
 *
 */
public class SQLiteConstruct extends SQLiteParent{

	protected SQLiteConstruct(Connection conn, int limitQueries) {
		super(conn, limitQueries);
	}

	protected SQLiteConstruct(Connection connection) {
		super(connection);
	}

	/**
	 * return a SQLiteConstruct object from a sqlite database
	 * with the connection to the database built in
	 * DON'T FORGET TO CALL close() METHOD WHEN YOU FINISHED TO USE SQLiteAccess
	 * @param fullPath - the path oh the database
	 * @return - SQLiteConstruct
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLiteAccess
	 */
	public static SQLiteConstruct getConnectionWithDatabase(String fullPath)
	throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException{
		Class.forName("org.sqlite.JDBC").newInstance();
		Connection conn = DriverManager.getConnection("jdbc:sqlite:/"+fullPath);
		conn.setAutoCommit(false);
		return new SQLiteConstruct(conn,100000);
	}
	/**
	 * return a SQLiteConstruct object from a sqlite database
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
	 * @throws SQLiteConstruct
	 */
	protected static SQLiteConstruct getConnectionWithDatabase(String fullPath,int limitQueriesSize)
	throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException{
		Class.forName("org.sqlite.JDBC").newInstance();
		Connection conn = DriverManager.getConnection("jdbc:sqlite:/"+fullPath);
		conn.setAutoCommit(false);
		return new SQLiteConstruct(conn,limitQueriesSize);
	}

	/**
	 * create a new database
	 * @param ext : must be qualitative or quantitative
	 */
	public void createNewDatabase(SQLiteExtension ext) throws SQLException{
		Statement stat;
		stat = this.connection.createStatement();
		stat.executeUpdate("drop table if exists attributes;");
		stat.executeUpdate("create table attributes (key text, value text);");
		PreparedStatement prep = this.connection.prepareStatement("insert into attributes values (?,?);");
		prep.setString(1, "datatype");
		prep.setString(2,ext.name());
		prep.executeUpdate();
		this.connection.commit();
	}

	/**
	 * create a new database
	 * @param type : must be qualitative or quantitative
	 */
	public void addAtribute(String key,String value) throws SQLException{
		Statement stat;
		stat = this.connection.createStatement();
		PreparedStatement prep = this.connection.prepareStatement("insert into attributes values (?,?);");
		prep.setString(1,key);
		prep.setString(2,value);
		prep.executeUpdate();
		this.connection.commit();
	}

	/**
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
	 * create a new table for a chromosome for 
	 * QUALITATIVE EXTENDED data
	 * initialize a prepared statement which you can insert values in the table
	 * (start,end,score,name,strand,attributes,type,id)
	 * @param chromosome - the chromosome name
	 * @throws SQLException
	 */
	public void newChromosome_qual_extended(String chromosome) throws SQLException {
		PreparedStatement stat = this.connection.prepareStatement("create table if not exists " +
				""+protect(chromosome)+" (start integer,end integer,score real,name text,strand integer,attributes text,type integer,id text);");
		stat.execute();
		this.connection.commit();
		this.insertStatement.put(
				chromosome,this.connection.prepareStatement("insert into "+protect(chromosome)+" values (?,?,?,?,?,?,?,?);"));

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
		System.out.println("get chrnames");
		List<String> chrNames = new ArrayList<String>();
		Statement stat = connection.createStatement();
		String query = "SELECT name FROM sqlite_master where type='table'and name!='attributes' and name!='types';";
		ResultSet rs = getResultSet(stat, query);
		while (rs.next()) {
			String name = rs.getString(1);
			System.out.println(name);
			chrNames.add(name);
		}
		rs.close();
		return chrNames;
	}
	/**
	 * look if the chromosome is already created
	 * @param chromosome - the chromosome name
	 * @return - true if the chromosome is created
	 */
	public boolean isCromosomeCreated(String chromosome){
		return insertStatement.get(chromosome)!=null;
	}


	/**
	 * write values in the insert statement of the qualitative database
	 * YOU MUST CALL newChromosome_qual(String chromosome) before
	 * @param chr - the chromosome
	 * @param start
	 * @param end
	 * @param score
	 * @param name - name of the feature
	 * @param strand
	 * @throws SQLException
	 */
	public void writeValues_qual(String chr, Integer start, Integer end,Float score,
			String name, Integer strand,String attributes) throws SQLException{
		PreparedStatement prep = insertStatement.get(chr);
		if(null==score){
			score=0f;
		}
		if(null==strand){
			strand=1;
		}
		if(null!=prep){
			prep.setInt(1,start);
			prep.setInt(2,end);
			prep.setFloat(3,score);
			prep.setString(4, name);
			prep.setInt(5, strand);
			prep.setString(6, attributes);
			prep.execute();
			verifyNbQueries();
		} else {
			System.err.println(
					"the prepared statement for the chromosome "+chr+" is null, you must call newChromosome_qual("+chr+") before ");
		}
	}
	/**
	 * write values in the insert statement of the qualitative database
	 * YOU MUST CALL newChromosome_qual(String chromosome) before
	 * @param chr - the chromosome
	 * @param start
	 * @param end
	 * @param score
	 * @param name - name of the feature
	 * @param strand
	 * @throws SQLException
	 */
	public void writeValues_qual_extended(String chr, Integer start, Integer end,Float score,
			String name, Integer strand,String attributes,Integer type,String id) throws SQLException{
		PreparedStatement prep = insertStatement.get(chr);
		if(null!=prep){
			prep.setInt(1,start);
			prep.setInt(2,end);
			prep.setFloat(3,score);
			prep.setString(4, name);
			prep.setInt(5, strand);
			prep.setString(6, attributes);
			prep.setInt(7, type);
			prep.setString(8, id);
			prep.execute();
			verifyNbQueries();
		} else {
			System.err.println(
					"the prepared statement for the chromosome "+chr+" is null, you must call newChromosome_qual_extended("+chr+") before ");
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
	public void writeValues_quant(String chr,Integer start,Integer end,Float score) throws SQLException{
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
	public void finalizeDatabase(Map<String,Integer> chr_length,boolean indexOnName,boolean indexOnScore,boolean indexOnRange) throws SQLException{
		Statement stat = this.connection.createStatement();
		stat.executeUpdate("create table chrNames (name text, length integer);");
		PreparedStatement prep = this.connection.prepareStatement("insert into chrNames values (?,?);");
		for(Map.Entry<String, Integer>entry : chr_length.entrySet()){
			String chr = entry.getKey();
			Statement stat2 = this.connection.createStatement();
			if(indexOnRange){
				stat2.executeUpdate("create index "+protect(chr+"_range_idx")+" on "+protect(chr)+" (start,end);");
			}
			if(indexOnName){
				stat2.executeUpdate("create index "+protect(chr+"_name_idx")+" on "+protect(chr)+" (name);");
			}
			if(indexOnScore){
				stat2.executeUpdate("create index "+protect(chr+"_score_idx")+" on "+protect(chr)+" (score);");
			}
			int length = entry.getValue();
			prep.setString(1, chr);
			prep.setInt(2,length);
			prep.executeUpdate();
		}
		this.connection.commit();
	}

	public void finalizeExtended_qual(List<String> types) throws SQLException {
		Statement stat = this.connection.createStatement();
		stat.executeUpdate("create table types (identifier integer, type text);");
		PreparedStatement prep = this.connection.prepareStatement("insert into types values (?,?);");
		for(int i=0;i<types.size();i++){
			prep.setInt(1, i);
			prep.setString(2, types.get(i));
			prep.execute();
		}

	}


	//FOR THE DAEMON <tranform_to_sqlite>
	/**
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
