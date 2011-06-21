package Peppy;

public class Property {
	
	private String name;
	private String description;
	private String value;
	
	private static int typeValue = 0;
	private static final int BOOLEAN = typeValue++;
	private static final int INT = typeValue++;
	private static final int DOUBLE = typeValue++;
	private static final int STRING = typeValue++;
	
	public Property(String name, String description) {
		this.name = name;
		this.description = description;
	}
	
	public String getName() {
		return name;
	}
	
	public String getDescription() {
		return description;
	}
	
	public boolean getBoolean() {
		return Boolean.parseBoolean(value);
	}
	
	public int getInt() {
		return Integer.parseInt(value);
	}
	
	public double getDouble() {
		return Double.parseDouble(value);
	}
	
	public String getString() {
		return value;
	}

}
