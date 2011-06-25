package ch.epfl.bbcf.bbcfutils.conversion.json.pojo;

import java.util.ArrayList;
import java.util.List;

import ch.epfl.bbcf.bbcfutils.conversion.json.Constants;

public class TrackDataQuantitative {

	private List<ZoomLevel> zoomLevels;
	private Integer tileWidth,min,max;
	
	
	public TrackDataQuantitative(int tileWidth,int min,int max,String databaseName,String chromosome){
		this.max=max;
		this.min=min;
		this.tileWidth=tileWidth;
		this.zoomLevels=new ArrayList<TrackDataQuantitative.ZoomLevel>();
		for(int zoom : Constants.zooms){
			ZoomLevel zl = new ZoomLevel(
					databaseName+"/"+chromosome+"_"+zoom+".db",
					Constants.JSON_HEIGHT,
					Integer.toString(zoom*100)
					);
			zoomLevels.add(zl);
		}
	}
	
	
	public static class ZoomLevel{
		public ZoomLevel(String url,String h,String bpt){
			this.urlPrefix=url;
			this.height=h;
			this.basesPerTile=bpt;
		}
		private String urlPrefix,height,basesPerTile;

		public void setUrlPrefix(String urlPrefix) {
			this.urlPrefix = urlPrefix;
		}

		public String getUrlPrefix() {
			return urlPrefix;
		}

		public void setHeight(String height) {
			this.height = height;
		}

		public String getHeight() {
			return height;
		}

		public void setBasePerTile(String basePerTile) {
			this.basesPerTile = basePerTile;
		}

		public String getBasePerTile() {
			return basesPerTile;
		}
	}






	public void setZoomLevels(List<ZoomLevel> zoomLevels) {
		this.zoomLevels = zoomLevels;
	}






	public List<ZoomLevel> getZoomLevels() {
		return zoomLevels;
	}









	public void setMax(Integer max) {
		this.max = max;
	}






	public Integer getMax() {
		return max;
	}






	public void setMin(Integer min) {
		this.min = min;
	}






	public Integer getMin() {
		return min;
	}






	public void setTileWidth(Integer tileWidth) {
		this.tileWidth = tileWidth;
	}






	public Integer getTileWidth() {
		return tileWidth;
	}
}
