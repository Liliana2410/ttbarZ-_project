from ROOT import TH1F


def histos_b():
	"""
	TH1F objetcs creation.
	- Convention of naming is: plot_kinematicOrTopologicalVariable_particle.
	- b1 an b2 are the leading b-jets and so on.
	Return: List of TH1F objetcs. 
	The last histo is plot_MET and the previous one is plot_PT_jets.
	"""
	#b-quark plots
	plot_PT_b1 = TH1F("PT_b1", "PT_b1", 100, 0.0, 1500.0)
	plot_PT_b2 = TH1F("PT_b2", "PT_b2", 100, 0.0, 1500.0)
	plot_P_b1 = TH1F("P_b1", "P_b1", 100, 0.0, 2000.0)
	plot_P_b2 = TH1F("P_b2", "P_b2", 100, 0.0, 2000.0)
	plot_PTvsPL = TH1F("PTvsPL", "PTvsPL", 100, 0.0, 3000.0)


	return (plot_PT_b1, plot_PT_b2, plot_P_b1, plot_P_b2, plot_PTvsPL)


