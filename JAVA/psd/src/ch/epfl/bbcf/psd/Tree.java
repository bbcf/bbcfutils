package ch.epfl.bbcf.psd;

import java.io.File;
import java.sql.ResultSet;
import java.sql.SQLException;

import org.apache.log4j.Logger;



public class Tree {

	private static final int TAB_WIDTH = Main.TAB_WIDTH;
	private Node leaf;
	private int imageNumber;
	private ConnectionStore connections;
	private String chromosome;
	private String out;

	private static Logger logger = Main.logger;

	public Tree(ConnectionStore connectionStore, String chromosome, String out) {
		logger.debug("tree initialization");

		this.out = out;
		this.connections = connectionStore;
		this.chromosome = chromosome;

		Node lastZoom = new Node(this, 100000, null, null, null, true);
		for(int i=4;i>=0;i--){
			Node node5 = new Node(this, (int)Math.pow(10, i) * 5, lastZoom, null, null, false);
			Node node2 = new Node(this, (int)Math.pow(10, i) * 2, null, null, null, false);
			lastZoom = new Node(this, (int)Math.pow(10, i), null, node2, node5, false);
		}
		this.leaf = lastZoom;
	}


	public void process(ResultSet scores) throws SQLException {
		logger.debug("processing scores");
		boolean hasScore = false;
		while(scores.next()){
			int start = scores.getInt(1);
			int stop = scores.getInt(2);
			float score = scores.getFloat(3);
			for(int j=start;j<=stop;j++){
				this.imageNumber = getImageNumber(j);
				leaf.fill(j, score, imageNumber);
			}
			hasScore = true;
		}

		if(hasScore){
			logger.debug("end writing");
			leaf.endWriting();
		} else {
			logger.debug("no score on chromosome " + this.chromosome +". Removing databases.");
			this.removeChromosome();
		}
	}

	private void removeChromosome() {
		File dir = new File(out);
		if (!dir.exists()) {
			logger.error(dir.getAbsolutePath() + " doesn't exist.");
			return;
		}

		String[] info = dir.list();

		for (int i = 0; i < info.length; i++) {
			File n = new File(dir + File.separator + info[i]);
			if (n.isFile() && info[i].startsWith(chromosome)) { 
				logger.debug("removing : "+ n.getPath());
				if(!n.delete()){
					logger.error("cannot delete "+n.getAbsolutePath());
				}
			}
		}
	}


	private static int getImageNumber(int position){
		double nb = (double)position / TAB_WIDTH;
		return (int)Math.ceil(nb);
	}


	public void writeValues(Node node) throws SQLException {
		this.writeValues(node, false);
	}

	public void writeValues(Node node, boolean finish) throws SQLException {
		SQLiteConnector.filldb(this.chromosome, node.getTab(), 
				node.getImageNumber(), node.getZoom(), this.connections, finish);

	}

}



