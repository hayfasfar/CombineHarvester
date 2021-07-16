import json

categories = [
    "mumu_OS", "mumu_SS",
    "ee_OS", "ee_SS",
    "mue_OS", "mue_SS",
    "emu_OS", "emu_SS",
]

category_dict = {
    1: "merged, prompt",
    2: "merged, intermediate",
    3: "merged, displaced",
    4: "resolved, prompt",
    5: "resolved, intermediate",
    6: "resolved, displaced",

}

rename_dict = {
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
    "lumi_2016": "luminosity (2016)",
    "lumi_2017": "luminosity (2017)",
    "lumi_2018": "luminosity (2018)"
}

nbins = 6
for category in categories:
    for region in ["A", "B", "C"]:
        for year in ["2016", "2017", "2018"]:
            for i in range(1, nbins+1):
                rename_dict["rate_bkg_{}_{}_bin{}_{}".format(category, region, i, year)] = "yield {}, region {}, {} bin, {}".format(category, region, category_dict[i], year)

with open("rename.json", "w+") as json_file:
    json.dump(rename_dict, json_file)


