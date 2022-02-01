import ROOT
from ROOT import TLorentzVector, TH1F
import histos_b
from collections import Counter
from math import sin, cos, sinh

def jet_reconstruction(jets, bs):

	global dijet_1
	global b_dijet_1
	global reconstructed_W1
	#merged category of the first pair of hadronic jets coming from the top.
	global notMerged_1
	global partiallyMerged_1
	global fullyMerged_1
	#index of the b-jet from the hadronic decaying top
	global index_b
	#index of the two jets coming from the hadronic top.
	global index_j1
	global index_j2
	#index of the two b-jets comming from the Z' decay.
	global index_bi
	global index_bj

	global Delta_PT_b_alpha_b_beta

	best_Err = 999999.9
	MW = 80.379
	Mt = 172.76
	smallest_pt_bibj = 999999.0
	
	for l1 in range(len(bs)):
		for l2 in range(l1 + 1,len(bs)):

			diff = abs(bs[l1][0].Pt() - bs[l2][0].Pt())
			if(diff < smallest_pt_bibj):

				smallest_pt_bibj = diff
				index_bi = l1
				index_bj = l2

	index_bi = 0
	index_bj = 1

	Delta_PT_b_alpha_b_beta = abs(smallest_pt_bibj)


	index_j_temp1 = -1
	index_j_temp2 = -1
	
	index_b_temp1 = -1
	
	
	#loop over the first group of particles (j_1, j_2, b_1).
	for i in range(len(jets)):
		for j in range(i + 1, len(jets)):
			for k in range(len(bs)):
				if ( (k != index_bi) and (k != index_bj) ):
					jtemp1 = jets[i]
					jtemp2 = jets[j]
					btemp1 = bs[k][0]
						
					Err = (abs((jtemp1 + jtemp2 + btemp1).M() - Mt))*100/Mt

					#Selection criteria.
					if (Err < best_Err):
						best_Err = Err
						index_j_temp1 = i
						index_j_temp2 = j
						index_b_temp1 = k
					
	
	#Initialize the first group of particles.				
	index_j1 = index_j_temp1
	index_j2 = index_j_temp2	
	index_b1 = index_b_temp1
	
	dijet_1 = jets[index_j1] + jets[index_j2]		
	dr_dijet_1 = jets[index_j1].DeltaR(jets[index_j2])
	b_dijet_1 = jets[index_j1] + jets[index_j2] + bs[index_b1][0]
	dr_b_dijet_1 = dijet_1.DeltaR(bs[index_b1][0])
    		
	#Merged category for the first group of particles.
	if (dr_dijet_1 > 0.8):

	 	notMerged_1 = True
	 	partiallyMerged_1 = False
	 	fullyMerged_1 = False


   	else:

		notMerged_1 = False
		
        	if (dr_b_dijet_1 > 1.0):
        
			partiallyMerged_1 = True
			fullyMerged_1 = False
			reconstructed_W1 = dijet_1

		else:

			partiallyMerged_1 = False
			fullyMerged_1 = True
		

	return jets[index_j1], jets[index_j2], bs[index_b1][0], bs[index_bi][0], bs[index_bj][0]


def frac_twobody_pTvspL(particle1, particle2):

	px1 = (particle1.Pt())*(cos(particle1.Phi()))
	py1 = (particle1.Pt())*(sin(particle1.Phi()))
	pz1 = (particle1.Pt())*(sinh(particle1.Eta()))
	
	px2 = (particle2.Pt())*(cos(particle2.Phi()))
	py2 = (particle2.Pt())*(sin(particle2.Phi()))
	pz2 = (particle2.Pt())*(sinh(particle2.Eta()))

	pT = ((px1 + px2)**2 + (py1 + py2)**2)**0.5
	pL = abs(pz1 + pz2)
	p = ((px1 - px2)**2 + (py1 - py2)**2 + (pz1 - pz2)**2)**0.5
	
	frac = pT/p

	return frac


def lepton_PT(TLV):
	return TLV[0].Pt()

def lepton_P(TLV):
	return TLV[0].P()

def PT(TLV):
	return TLV.Pt()

def P(TLV):
	return TLV.P()
    
def histos_fill(plot, variable):
	return plot.Fill(variable)


def histos_Draw(plot):
	return plot.Draw('HISTOS')  


def histos_Write(plot):
	return plot.Write()


def histos_Reset(plot):
	return plot.Reset('ICESM') 


#signals = ["ttbarh", "ttbarbbar_noh", "ttbarttbar", "Zprime_bbar_350", "Zprime_bbar_3000"]
#jobs = [2,2,2,2,2]

signals = ["ttbarh", "ttbarbbar_noh", "ttbarttbar", "Zprime_bbar_3000", "Zprime_bbar_1000", "Zprime_bbar_350"]
jobs = [1,1,1,1,1,1]

#------------------ HISTOGRAMS ---------------------
c1 = ROOT.TCanvas("c1", "Titulo")    # ROOT canvas

plots = histos_b.histos_b()   # plots is a list of TH1F objects

plot_P_genb1 = TH1F("P_genb1", "P_genb1", 100, 0.0, 2500.0)
plot_P_genb2 = TH1F("P_genb2", "P_genb2", 100, 0.0, 2500.0)
plot_PT_genb1 = TH1F("PT_genb1", "PT_genb1", 100, 0.0, 2000.0)
plot_PT_genb2 = TH1F("PT_genb2", "PT_genb2", 100, 0.0, 2000.0)
plot_ETA_genb1 = TH1F("ETA_genb1", "ETA_genb1", 100, -5, 5)
plot_ETA_genb2 = TH1F("ETA_genb2", "ETA_genb2", 100, -5, 5)
plot_PHI_genb1 = TH1F("PHI_genb1", "PHI_genb1", 100, -4, 4)
plot_PHI_genb2 = TH1F("PHI_genb2", "PHI_genb2", 100, -4, 4)

plot_DeltaR_genb1b2 = TH1F("DeltaR_genb1b2", "DeltaR_genb1b2", 100, -4, 4)
plot_DeltaPhi_genb1b2 = TH1F("DeltaPhi_genb1b2", "DeltaPhi_genb1b2", 100, -4, 4)
plot_gen_fromzp = TH1F("gen_fromzp", "gen_fromzp", 2, 0.0, 3.0)

#------------- ITERATING THE FILES AND MAKING THE HISTOGRAMS ----------------

for n_signal, signal in enumerate(signals):

	#List of PDG-id of first mother of gen bs that have a m1 in all events.
	b_mother1 = []
	#List of PDG-id of first daughter of gen bs that have a d1 in all events.
	b_daughter1 = []
	
	#Number of b-quarks with a daughter in all events.
	withdaughter = 0
	#Number of b-quarks with a mother in all events.
	withmother = 0
	#Number of charms with a mother in all events.
	charmwithmother = 0
	#Number of charms with a b mother in all events.
	fromb = 0
	#Number of b-quark with a Z' as its mother in all events.
	fromzp = 0
	#Number of b-quarks coming from a Z' with daughter.
	fromzpwithdaughter = 0
	#Number of b-quarks coming from a Z' that also have a neutrino as a daughter.
	toneutrino = 0
	#Counters of times the chosen bs from the Z' (bi - bj) match with gen bs coming from a Z' (smallest DR) in all events.
	match1 = 0
	match2 = 0
	
	#Number of good events (match1 and match2 >= 1) in the non-problematic region of b2: pt(b2) > 100GeV.
	right_events = 0
	#Number of good events (match1 and match2 >= 1) in the problematic region of b2: 20GeV < pt(b2) < 100GeV.
	right_events_p = 0
	#Number of events that pass the cross-cleaning in the non-problematic region of b2.
	total_events = 0
	#Number of events that pass the cross-cleaning in the problematic region of b2.
	total_events_p = 0

	f = ROOT.TFile(signal + ".root", "recreate")

	for ind in range(1, jobs[n_signal] + 1):

		directory = str("/disco4/SIMULACIONES/from_disco3/with_delphes/" + signal + "/" + signal + "_" + str(ind) + "/Events/run_01/tag_1_delphes_events.root")
		File = ROOT.TChain("Delphes;1")
		File.Add(directory)
		Number = File.GetEntries()

		print("Signal: " + signal + "_" + str(ind))

		for i in range(Number):
			Entry = File.GetEntry(i)

			#Initializes particles lists.
			jets = []
			bs = []
			METs = []
			taus = []
			electrons = []
			muons = []
			

			EntryFromBranch_j = File.Jet.GetEntries()
			for j in range(EntryFromBranch_j):

				BTag = File.GetLeaf("Jet.BTag").GetValue(j)
				TauTag = File.GetLeaf("Jet.TauTag").GetValue(j)

				#searches for jets.
				if (BTag != 1 and TauTag != 1):
					jet = TLorentzVector()
					jet_PT, jet_Eta, jet_Phi, jet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
					jet.SetPtEtaPhiM(jet_PT, jet_Eta, jet_Phi, jet_M)
					jets.append(jet)

				#searches for b_jets.
				elif (BTag == 1 and TauTag != 1):
					bjet = TLorentzVector()
					bjet_PT, bjet_Eta, bjet_Phi, bjet_M, bjet_D0 = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j), File.GetLeaf("Jet.Charge").GetValue(j)
					bjet.SetPtEtaPhiM(bjet_PT, bjet_Eta, bjet_Phi, bjet_M)
					#Creates the list of bs as a touple containing the normal TLorentzVector and their impact parameters as it is needed.
					bs.append((bjet, bjet_D0))

				#searches for taus.
				elif (TauTag == 1 and BTag != 1):
					tau = TLorentzVector()
					tau_PT, tau_Eta, tau_Phi, tau_M, tau_Charge = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j),  File.GetLeaf("Jet.Charge").GetValue(j)
					tau.SetPtEtaPhiM(tau_PT, tau_Eta, tau_Phi, tau_M)
					taus.append(tau)


			# MET (neutrinos).
			Total_MET = 0
			EntryFromBranch_MET = File.MissingET.GetEntries()
			for j in range(EntryFromBranch_MET):
				MET = TLorentzVector()
				MET_PT, MET_Eta, MET_Phi, MET_M  = File.GetLeaf("MissingET.MET").GetValue(j), File.GetLeaf("MissingET.Eta").GetValue(j), File.GetLeaf("MissingET.Phi").GetValue(j), 0.0
				MET.SetPtEtaPhiM(MET_PT, MET_Eta, MET_Phi, MET_M)
				METs.append(MET)
				Total_MET += MET_PT

			EntryFromBranch_e = File.Electron.GetEntries()
			for j in range(EntryFromBranch_e):
				electron = TLorentzVector()
				electron_PT, electron_Eta, electron_Phi, electron_M  = File.GetLeaf("Electron.PT").GetValue(j), File.GetLeaf("Electron.Eta").GetValue(j), File.GetLeaf("Electron.Phi").GetValue(j), 520.998
				electron.SetPtEtaPhiM(electron_PT, electron_Eta, electron_Phi, electron_M)
				electrons.append(electron)

			EntryFromBranch_mu = File.Muon.GetEntries()
			for j in range(EntryFromBranch_mu):
				muon = TLorentzVector()
				muon_PT, muon_Eta, muon_Phi, muon_M  = File.GetLeaf("Muon.PT").GetValue(j), File.GetLeaf("Muon.Eta").GetValue(j), File.GetLeaf("Muon.Phi").GetValue(j), 105658.374
				muon.SetPtEtaPhiM(muon_PT, muon_Eta, muon_Phi, muon_M)
				muons.append(muon)

			leptons = electrons + muons
      			#Checks if the event has the minimum theoretical expected number of particles (the count might be higher due to particle-detector effects).
			if (len(jets) >= 2 and len(bs) >= 4 and len(leptons) != 0):


				jets.sort(reverse = True, key=PT)     
				bs.sort(reverse = True, key=lepton_P)
				taus.sort(reverse = True, key=PT)


				#Defines the real jets and bs from the event.
				j1, j2, b, b1, b2 = jet_reconstruction(jets, bs)
				realjets = [j1, j2]
				realbs = [b, b1, b2]
				
				#List of all gen bs coming from a zp.
				gen_bs = []
				#List of all gen neutrinos that come from a gen b that is the daughter of a zp.
				gen_neutrinos = []
				#List of all gen bs.
				gen_allbs = []
				#Indexes of the position of genb1 and genb2 in the array of genParticles.
				index_genb1 = 0
				index_genb2 = 0
				
				#Counters of times the chosen bs (bi - bj) match with gen-bs (smallest DR) from a Z' in the event (if it is in the non-problematic region).
				temp_match1 = 0
				temp_match2 = 0

				#Counters of times the chosen bs (bi - bj) match with gen-bs (smallest DR) from a Z' in the event (if it is in the problematic region).
				temp_match1_p = 0
				temp_match2_p = 0
				
				
				#~CODE FOR GEN-BS AND MATCHING ALGORITHM.
				EntryFromBranch_gen = File.Particle.GetEntries()
				for j in range(EntryFromBranch_gen):

					gen_id = File.GetLeaf("Particle.PID").GetValue(j)

					#Checks if the GenParticle is a charm.
					if (abs(gen_id) == 4):
						#Checks amount of charms with a b mother.
						if (File.GetLeaf("Particle.M1").GetValue(j) >= 0.0):
							charmwithmother += 1
							if ( abs(File.GetLeaf("Particle.PID").GetValue(int(File.GetLeaf("Particle.M1").GetValue(j)))) == 5):
								fromb += 1


					#Checks if the GenParticle is a b.
					if (abs(gen_id) == 5):
					
						genb = TLorentzVector()
						genb_PT, genb_Eta, genb_Phi, genb_M  = File.GetLeaf("Particle.PT").GetValue(j), File.GetLeaf("Particle.Eta").GetValue(j), File.GetLeaf("Particle.Phi").GetValue(j), File.GetLeaf("Particle.Mass").GetValue(j)
						genb.SetPtEtaPhiM(genb_PT, genb_Eta, genb_Phi, genb_M)


						#Organizes the array gen_allbs, array of all gen-bs, by P, while keeping track of the index in the genFile of the two most energetic bs.
						if (len(gen_allbs) == 0):
							gen_allbs.append(genb)
							index_genb1 = j

						elif (len(gen_allbs) == 1):
							if (genb.P() >= gen_allbs[0].P()):
								gen_allbs.insert(0, genb)
								index_genb1 = j
							else:
								gen_allbs.insert(1, genb)
								index_genb2 = j

						else:
							added = False
							l = 0
							while (added == False):
								if (l < len(gen_allbs) - 1):

									if (genb.P() >= gen_allbs[0].P()):
										gen_allbs.insert(0, genb)
										index_genb1 = j
										added = True

									elif (genb.P() < gen_allbs[l].P()) and (genb.P() >= gen_allbs[l + 1].P()):
										gen_allbs.insert(l + 1, genb)
										added = True
										if (l == 0):
											index_genb2 = j

								else:
									if (genb.P() >= gen_allbs[-1].P()):
										gen_allbs.insert(-1, genb)
										added = True
									else:
										gen_allbs.append(genb)
										added = True

								l += 1
								

						#Checks if the tau is in the problematic region of taus.
						#if ((genjet.Pt() > 15) and (genjet.Pt() < 100)):
							#if (abs(genjet.Eta()) < 2.3):
						
						#Mother and grandmother index in GenParticle array, < 0 if there is no such particle.
						mother1 = File.GetLeaf("Particle.M1").GetValue(j)
						mother2 = File.GetLeaf("Particle.M2").GetValue(j)
						
						if (mother1 >= 0.0):

							withmother += 1
							#Mother ID of the gen-b.
							mother1_id =  File.GetLeaf("Particle.PID").GetValue(int(mother1))
							b_mother1.append(abs(mother1_id))

							#Checks if the b has the zp as its mother.
							if (abs(mother1_id) == 32):

								fromzp += 1
								
								daughter1 = File.GetLeaf("Particle.D1").GetValue(j)
								if (daughter1 >= 0.0):
									
									fromzpwithdaughter += 1
									daughter1_id =  File.GetLeaf("Particle.PID").GetValue(int(daughter1))
									b_daughter1.append(abs(daughter1_id))	
									
									#Checks if the b has a neutrino as its daughter.
									if ((abs(daughter1_id) == 12) or (abs(daughter1_id) == 14) or (abs(daughter1_id) == 16)):
										
										toneutrino += 1	
										neutrino = TLorentzVector()
										neutrino_PT, neutrino_Eta, neutrino_Phi, neutrino_M  = File.GetLeaf("Particle.PT").GetValue(int(daughter1)), File.GetLeaf("Particle.Eta").GetValue(int(daughter1)), File.GetLeaf("Particle.Phi").GetValue(int(daughter1)), 0
										neutrino.SetPtEtaPhiM(neutrino_PT, neutrino_Eta, neutrino_Phi, neutrino_M)
										gen_neutrinos.append(neutrino)
					
								gen_bs.append(genb)
								
								#Index and dR of the detector-b that is closest to the gen-b coming from a Z'.
								best_b_index = -1
								smallest_DR= 999999.9

								for i in range (len(bs)):
									DR = genb.DeltaR(bs[i][0])
									if (DR < smallest_DR):
										smallest_DR = DR
										best_b_index = i

								
								#Checks if the event is in the non-problematic region.
								if (b2.Pt() >= 100):
									#Checks if the pion matched with tau1.
									#if (best_b_index == index_bi):
									if (genb.DeltaR(b1) <= 0.3):
										match1 += 1
										temp_match1 += 1

									#Checks if the pion matched with tau2.
									#elif (best_b_index == index_bj):
									elif (genb.DeltaR(b2) <= 0.3):
										match2 += 1
										temp_match2 += 1

								#Checks if the event is in the problematic region.
								if (b2.Pt() < 100):
									#Checks if the pion matched with tau1.
									#if (best_b_index == index_bi):
									if (genb.DeltaR(b1) <= 0.3):
										temp_match1_p += 1
									#Checks if the pion matched with tau2.
									#elif (best_b_index == index_bj):
									elif (genb.DeltaR(b2) <= 0.3):
										temp_match2_p += 1
		
		
				
					
				#Checks if the event had the minimun required gen-pions coming from a gen-tau as to count as possible correctly tagged event.
				if (len(gen_bs) >= 2):
				
					#Checks if both taus from the Z' matched with a gen-pion in the non-problematic region.
					if ((temp_match1 >= 1) and (temp_match2 >= 1)):
						right_events += 1

					#Checks if both taus from the Z' matched with a gen-pion in the problematic region.
					elif ((temp_match1_p >= 1) and (temp_match2_p >= 1)):
						right_events_p += 1
						
					#Checks if the event was in the non-problematic region.
					if (b2.Pt() >= 100):
						total_events += 1

					#Checks if the event was in the problematic region.
					elif (b2.Pt() < 100):
						total_events_p += 1	
							

				#gen_allbs.sort(reverse = True, key=P)

				genb1 = gen_allbs[0]
				genb2 = gen_allbs[1]

				genb1mother = int(File.GetLeaf("Particle.M1").GetValue(index_genb1))
				genb2mother = int(File.GetLeaf("Particle.M1").GetValue(index_genb2))

				histos_fill(plot_P_genb1, genb1.P())
				histos_fill(plot_P_genb2, genb2.P())
				histos_fill(plot_PT_genb1, genb1.Pt())
				histos_fill(plot_PT_genb2, genb2.Pt())
				histos_fill(plot_ETA_genb1, genb1.Eta())
				histos_fill(plot_ETA_genb2, genb2.Eta())
				histos_fill(plot_PHI_genb1, genb1.Phi())
				histos_fill(plot_PHI_genb2, genb2.Phi())

				histos_fill(plot_DeltaR_genb1b2, genb1.DeltaR(genb2))
				histos_fill(plot_DeltaPhi_genb1b2, genb1.DeltaPhi(genb2))

				if (genb1mother >= 0.0) and (genb2mother >= 0.0):

					pdggenb1mother = File.GetLeaf("Particle.PID").GetValue(genb1mother)
					pdggenb2mother = File.GetLeaf("Particle.PID").GetValue(genb2mother)

					if ((pdggenb1mother == 32) and (pdggenb2mother == 32)):

						plot_gen_fromzp.AddBinContent(1)

					else:
						plot_gen_fromzp.AddBinContent(2)
			
				else:
					plot_gen_fromzp.AddBinContent(2)

				a = plot_gen_fromzp.GetXaxis()
				a.SetBinLabel(1,"Both from Z'")
				a.SetBinLabel(2,"Not both from Z'")

				#Pt requirement of the different jets (20 for taus is LHC minimum).
				if ((b1.Pt() > 30) and (b2.Pt() > 30)):

					#if (((temp_match1 >= 1) and (temp_match2 >= 1)) or ((temp_match1_p >= 1) and (temp_match2_p >= 1))):

					#Variables to be plotted.
					variables = [b1.Pt(), b2.Pt(), b1.P(), b2.P(), frac_twobody_pTvspL(b1, b2)]

					for i in range(len(plots)):
				  		histos_fill(plots[i], variables[i])


	c = Counter(b_mother1)
	c.most_common(1)
	print ("bs most common mother: ", len(b_mother1), c.most_common(1))

	c2 = Counter(b_daughter1)
	c2.most_common(1)
	print ("bs most common daughter: ", len(b_daughter1), c2.most_common(1))
	
	print ("charms with a mother: ", charmwithmother, "charms with a b mother: ", fromb, "bs with mother: ", withmother, "bs with Z' as a mother: ", fromzp, "bs coming from the Z' with a daughter: ", fromzpwithdaughter, "bs coming from the Z' with a neutrino as a daughter: ", toneutrino, "bi's that match with a b-quark: ", match1, "bj's that match with a b-quark: ", match2, "total events that pass cross-cleaning in the non-problematic region: ", total_events, "events with bi and bj matching in the non-problematic region: ", right_events, "total events that pass cross cleaning in the problematic region: ", total_events_p, "events with bi and bj matching in the problematic region: ", right_events_p)


	plot_P_genb1.Draw('HISTOS')
	plot_P_genb2.Draw('HISTOS')	
	plot_PT_genb1.Draw('HISTOS')
	plot_PT_genb2.Draw('HISTOS')
	plot_ETA_genb1.Draw('HISTOS')
	plot_ETA_genb2.Draw('HISTOS')
	plot_PHI_genb1.Draw('HISTOS')
	plot_PHI_genb2.Draw('HISTOS')

	plot_DeltaR_genb1b2.Draw('HISTOS')
	plot_DeltaPhi_genb1b2.Draw('HISTOS')
	plot_gen_fromzp.Draw('HISTOS')

	for plot in plots:
		histos_Draw(plot)

	# Updating the canvas
	c1.Update()

	plot_P_genb1.Write()
	plot_P_genb2.Write()
	plot_PT_genb1.Write()
	plot_PT_genb2.Write()
	plot_ETA_genb1.Write()
	plot_ETA_genb2.Write()
	plot_PHI_genb1.Write()
	plot_PHI_genb2.Write()

	plot_DeltaR_genb1b2.Write()
	plot_DeltaPhi_genb1b2.Write()
	plot_gen_fromzp.Write()

	for plot in plots:
		histos_Write(plot)

	# Closing the ROOT file where the histos were saved.
	f.Close()

	plot_P_genb1.Reset('ICESM')
	plot_P_genb2.Reset('ICESM')
	plot_PT_genb1.Reset('ICESM')
	plot_PT_genb2.Reset('ICESM')
	plot_ETA_genb1.Reset('ICESM')
	plot_ETA_genb2.Reset('ICESM')
	plot_PHI_genb1.Reset('ICESM')
	plot_PHI_genb2.Reset('ICESM')

	plot_DeltaR_genb1b2.Reset('ICESM')
	plot_DeltaPhi_genb1b2.Reset('ICESM')
	plot_gen_fromzp.Reset('ICESM')

	for plot in plots:
		histos_Reset(plot)

