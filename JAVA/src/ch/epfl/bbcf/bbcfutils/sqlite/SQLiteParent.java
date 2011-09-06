package ch.epfl.bbcf.bbcfutils.sqlite;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.Map;

public class SQLiteParent {

	/**
	 * a sqlite connection to a database
	 */
	protected Connection connection;
	/**
	 * key = chromosome name,so the table name
	 * value = a prepared statement
	 */
	protected Map<String,PreparedStatement> insertStatement;

	/**
	 * current number of queries
	 */
	protected int nbQueries;

	/**
	 * the number limit of queries to reach
	 * before committing 
	 */
	protected int limitQueries;


	protected SQLiteParent(Connection connection){
		this.connection = connection;
		this.insertStatement = new HashMap<String, PreparedStatement>();
		this.limitQueries = 100000;
		this.nbQueries = 0;
	}
	protected SQLiteParent(Connection conn,int limitQueries){
		this.connection = conn;
		this.insertStatement = new HashMap<String,PreparedStatement>();
		this.limitQueries = limitQueries;
		this.nbQueries=0;
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
	protected static ResultSet getResultSet(Statement statement,String query) throws SQLException{
		return statement.executeQuery(query);
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
	protected static String protect(String str){
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
	public void verifyNbQueries() throws SQLException{
		this.nbQueries++;
		if(nbQueries>limitQueries){
			this.commit();
			nbQueries = 0;
		}
	}

	
}
