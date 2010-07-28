package HMMScore;


public class IonMatch {
	private int ionName, ionIndex;
	double NTerCleavage, CTerCleavage;

	public IonMatch(int name, int index, double NCleavage, double CCleavage) {
		this.ionName = name;
		this.ionIndex = index;
		this.NTerCleavage = NCleavage;
		this.CTerCleavage = CCleavage;
	}

	public double getNTerCleavage() {
		return NTerCleavage;
	}

	public double getCTerCleavage() {
		return CTerCleavage;
	}

	public int getIonName() {
		return ionName;
	}

	public int getIonIndex() {
		return ionIndex;
	}

	public void setNTerCleavage(double _cleavage) {
		NTerCleavage = _cleavage;
	}

	public void setCTerCleavage(double _cleavage) {
		CTerCleavage = _cleavage;
	}

	public void setIonName(int _name) {
		ionName = _name;
	}

	public void setIonIndex(int _index) {
		ionIndex = _index;
	}

}
