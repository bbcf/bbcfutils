package ch.epfl.bbcf.psd;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ConnectionStore {

	private Map<String,Integer> indexes;
	private List<Connection> connections;
	private List<PreparedStatement> statements;
	private List<Integer> queries;
	
	public ConnectionStore() {
		indexes = new HashMap<String, Integer>();
		connections = new ArrayList<Connection>();
		statements = new ArrayList<PreparedStatement>();
		queries = new ArrayList<Integer>();
	}

	public void addDatabase(String database, Connection conn,
			PreparedStatement prep) {
		connections.add(conn);
		statements.add(prep);
		queries.add(0);
		indexes.put(database,queries.size()-1);
	}
	
	public Connection getConnection(String database){
		return connections.get(indexes.get(database));
	}
	public PreparedStatement getPreparedStatement(String database){
		return statements.get(indexes.get(database));
	}
	public void setPreparedStatement(String database, PreparedStatement prep) {
		statements.set(indexes.get(database), prep);
	}
	public int getNbQueries(String database){
		return queries.get(indexes.get(database));
	}

	public void setNbQueries(String database, int i) {
		queries.set(indexes.get(database), i);
	}

	public void destruct() throws SQLException {
		for(Connection c : connections){
			c.commit();
			c.close();
		}
		
	}

	
	
}
