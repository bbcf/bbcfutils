package ch.epfl.bbcf.bbcfutils.json;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.util.List;
import java.util.Map;

import org.codehaus.jackson.JsonProcessingException;
import org.codehaus.jackson.annotate.JsonAutoDetect;
import org.codehaus.jackson.map.DeserializationConfig;
import org.codehaus.jackson.map.ObjectMapper;
import org.codehaus.jackson.map.introspect.VisibilityChecker;
import org.codehaus.jackson.node.TreeTraversingParser;
import org.codehaus.jackson.type.TypeReference;

import ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo.GenrepObject;


public class JsonMapper
{
	private static final int INITIAL_SIZE = 2048;
	/** See http://wiki.fasterxml.com/JacksonBestPracticeThreadSafety?highlight=(\bCategoryJackson\b) */
	private static ObjectMapper mapper;

	static{
		mapper = new ObjectMapper();
		mapper.configure(DeserializationConfig.Feature.FAIL_ON_UNKNOWN_PROPERTIES,false);
		DeserializationConfig deserializationConfig = mapper.getDeserializationConfig();
		deserializationConfig.enable(DeserializationConfig.Feature.AUTO_DETECT_FIELDS);
		mapper.setVisibilityChecker(VisibilityChecker.Std.defaultInstance().withFieldVisibility(JsonAutoDetect.Visibility.ANY));
	}

	public static <T> String serialize(T o) throws IOException{
		StringWriter sw = new StringWriter(INITIAL_SIZE);
		mapper.writeValue(sw, o);
		return sw.toString();
	}

	public static <T extends GenrepObject> GenrepObject deserializeObject(String source, Class<? extends GenrepObject> targetClass) throws IOException{
		ByteArrayInputStream stream = new ByteArrayInputStream(source.getBytes());
		TreeTraversingParser treeTraversingParser = new TreeTraversingParser(mapper.readTree(stream));
		treeTraversingParser.setCodec(mapper);
		return treeTraversingParser.readValueAs(targetClass);
	}
	public static <K,V> Map<K,V> deserialize(String source, TypeReference<Map<K,V>> target) throws IOException{
		ByteArrayInputStream stream = new ByteArrayInputStream(source.getBytes());
		TreeTraversingParser treeTraversingParser = new TreeTraversingParser(mapper.readTree(stream));
		treeTraversingParser.setCodec(mapper);
		return treeTraversingParser.readValueAs(target);
	}

	public static <T extends GenrepObject> List<? extends GenrepObject> deserializeArray(String source,
			TypeReference<List<GenrepObject>> typeReference) throws JsonProcessingException, IOException {
		ByteArrayInputStream stream = new ByteArrayInputStream(source.getBytes());
		TreeTraversingParser treeTraversingParser = new TreeTraversingParser(mapper.readTree(stream));
		treeTraversingParser.setCodec(mapper);
		return treeTraversingParser.readValueAs(typeReference);
	}







}