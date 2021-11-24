import json

categories = [
    "mumu_OS", "mumu_SS",
    "ee_OS", "ee_SS",
    "mue_OS", "mue_SS",
    "emu_OS", "emu_SS",
]

def make_pretty(s):

    """ Make pretty the category text"""
    s = s.replace("mumu_OS", "\mu^{\pm}\mu^{\mp}")
    s = s.replace("mue_OS", "\mu^{\pm}e^{\mp}")
    s = s.replace("emu_OS", "e^{\pm}\mu^{\mp}")
    s = s.replace("ee_OS", "e^{\pm}e^{\mp}")
    s = s.replace("mumu_SS", "\mu^{\pm}\mu^{\pm}")
    s = s.replace("mue_SS", "\mu^{\pm}e^{\pm}")
    s = s.replace("emu_SS", "e^{\pm}\mu^{\pm}")
    s = s.replace("ee_SS", "e^{\pm}e^{\pm}")
    return s

category_dict = {
    1: "boosted, prompt",
    2: "boosted, intermediate",
    3: "boosted, displaced",
    4: "resolved, prompt",
    5: "resolved, intermediate",
    6: "resolved, displaced",

}

rename_dict_year = {
    "r" : "#sigma",
    "abcd_uncertainty" : "ABCD non-closure",
    "displaced_track": "displaced track",
    "pu" : "pileup",
    "jer": "jet energy resolution",
    "jesTotal": "jet energy scale",
    "unclEn": "unclustered energy",
    "loose_electron_reco": "displaced electron reconstruction",
    "tight_electron_reco": "prompt electron reconstruction",
    "tight_electron_id": "prompt electron identification",
    "tight_muon_id": "prompt muon identification",
    "tight_muon_iso": "prompt muon isolation",
    "trigger": "prompt muon trigger efficiency",
    "lumi_": "luminosity",
    "unc_boosted_": "non-closure boosted",
    "unc_resolved_": "non-closure resolved",
    "pdf": "pdf",
    "scale": "scale",
    "displaced_lepton_": "displaced lepton"
}

rename_dict = {}

for key, item in rename_dict_year.iteritems():
    for year in ["2016", "2017", "2018"]:
        rename_dict[key+year] = "{} ({})".format(item, year)

nbins = 6
for category in categories:
    for region in ["A", "B", "C", "D"]:
        for year in ["2016", "2017", "2018"]:
            for i in range(1, nbins+1):
                if region == "D":
                    rename_dict["rate_bkg_{}_{}_bin{}_{}_unc".format(category, region, i, year)] = "unc. {}, region {}, {} bin, {}".format(make_pretty(category), region, category_dict[i], year)
                else:
                    rename_dict["rate_bkg_{}_{}_bin{}_{}".format(category, region, i, year)] = "yield {}, region {}, {} bin, {}".format(make_pretty(category), region, category_dict[i], year)


with open("rename.json", "w+") as json_file:
    json.dump(rename_dict, json_file)


