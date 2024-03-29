
from cfg import *
sys.path.append('libraries')

from utils import *
from NeuronMorphology import NeuronMorphology
import convertAllenSpace as CAS


class NeuronPopulation:

    def __init__(self, data_path, res):

        annotation, neurite_length, acr_to_morpho_id, ancestorsById, id2acr = load_useful_variables(data_path, res)

        self.un_num = lambda x : x.split('1')[0].split('2/3')[0].split('4')[0].split('5')[0].split('6a')[0].split('6b')[0]
        self.res = res
        self.data_path = data_path
        self.annotation = annotation
        self.id2acr = id2acr
        self.acr2id = {val:key for key,val in id2acr.items()}
        self.acr_to_morpho_id = acr_to_morpho_id
        self.neurite_length = neurite_length
        self.ancestorsById = ancestorsById


    def has_numbers(self, inputString):
        return any(char.isdigit() for char in inputString)

    def make_connectivity_matrix(self, source_areas, target_areas, feature = 'counts', extract = 'terminals'):

        mesoscale_stats = {}
        targets_per_neuron = OrderedDict()
        somata  = []

        acr_to_morpho_id = deepcopy(self.acr_to_morpho_id)

        annot_shape = self.annotation.shape
        out_orientation = ['um({})'.format(self.res),'PIR','corner']
        if self.neurite_length is not None:
            nld_list = self.neurite_length.keys()
        for source_area in source_areas:
            if source_area not in acr_to_morpho_id:
                acr_to_morpho_id[source_area] = []
            for rand_area in acr_to_morpho_id.keys():  # First ensure that you have aggregated neurons of every child node
                if rand_area == '[background]' or rand_area == source_area: continue
                nu_id = str(self.acr2id[rand_area])
                if self.acr2id[source_area] in self.ancestorsById[nu_id]: # found a child area ...
                      acr_to_morpho_id[source_area].extend(acr_to_morpho_id[rand_area])
            if acr_to_morpho_id[source_area] == []: continue # nothing found here ...
            for neuron_id, neuron_soma_cord in acr_to_morpho_id[source_area]:
                # for each morphology, extract useful properties ...

                neuron_id = neuron_id.split('.')[0]
                nld_id = braintell_2_nld(nld_list, neuron_id)
                if nld_id == -1:
                    continue
                if self.un_num(self.neurite_length[nld_id]['soma']['region']) != self.un_num(self.neurite_length[nld_id]['soma']['correctedRegion']):
                    continue

                target_matches = [val for val in target_areas for val2 in self.neurite_length[nld_id]['axon'].keys() if val in val2]
                if len(target_matches) == 0: # does not target any of the areas
                    continue

                db = 'mouselight' if 'AA' in neuron_id else 'braintell'
                neuronPath = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}.json.gz'.format(db,neuron_id)
                file_content = requests.get(neuronPath).content
                target_neuron = json.loads(zlib.decompress(file_content, 16+zlib.MAX_WBITS))

                neuron_cls = NeuronMorphology(neuronDict = target_neuron)
                neuron_name = neuron_id.split('.')[0]
                if db == 'braintell':
                    re_str = neuron_name.split('_')[1:3]
                    neuron_name = re_str[0]+'_'+ re_str[1]
                neuron_cls.transform(out_orientation)

                anatomical_stats = neuron_cls.get_anatomical_stats(self.annotation,self.id2acr, mode = extract)

                if feature == 'length':
                    if extract == 'terminals':
                        term = True
                    else:
                        term = False
                    axonal_length = neuron_cls.compute_axonal_length(self.annotation,self.id2acr, terminals = term)
                    for key,val in axonal_length.items():
                        if key == '[background]' : continue
                        trg_id = str(self.acr2id[key])
                        trg_hits = [trg_val for trg_val in target_areas if self.acr2id[trg_val] in self.ancestorsById[trg_id]] # found a child area
                        if len(trg_hits) > 0:
                                key = trg_hits[0]
                                if (source_area,neuron_name) not in mesoscale_stats.keys():
                                    mesoscale_stats[(source_area,neuron_name)] = OrderedDict()
                                if key not in mesoscale_stats[(source_area,neuron_name)].keys():
                                    mesoscale_stats[(source_area,neuron_name)][key] = np.sum(val)
                                else:
                                    mesoscale_stats[(source_area,neuron_name)][key] += np.sum(val)
                elif feature == 'counts':
                    for key,val in anatomical_stats.items():
                        if key == '[background]' : continue
                        trg_id = str(self.acr2id[key])
                        trg_hits = [trg_val for trg_val in target_areas if self.acr2id[trg_val] in self.ancestorsById[trg_id]] # found a child area
                        if len(trg_hits) > 0:
                            key = trg_hits[0]
                            if (source_area,neuron_name) not in mesoscale_stats.keys():
                                mesoscale_stats[(source_area,neuron_name)] = OrderedDict()
                            if key not in mesoscale_stats[(source_area,neuron_name)].keys():
                                mesoscale_stats[(source_area,neuron_name)][key] = len(val)
                            else:
                                mesoscale_stats[(source_area,neuron_name)][key] += len(val)

                targets_per_neuron[neuron_name] = np.array([val for key,vals in anatomical_stats.items() for val in vals])
                soma_allen = np.round(neuron_cls.points[neuron_cls.somaPointIdx][0:3]).astype(int)
                somata.append(soma_allen)

        mesoscale_stats_df = pd.DataFrame(mesoscale_stats).T
        if len(list(mesoscale_stats.keys())) > 0:
            mesoscale_stats_df.index = mesoscale_stats_df.index.rename(['Source','Neuron Id'])
            re_sort_cols = np.argsort([col[0] for col in mesoscale_stats_df.columns])
            re_sort_rows = np.argsort([col[0] for col in mesoscale_stats_df.index])
            mesoscale_stats_df = mesoscale_stats_df.reindex(mesoscale_stats_df.columns[re_sort_cols], axis=1)
            mesoscale_stats_df = mesoscale_stats_df.reindex(mesoscale_stats_df.index[re_sort_rows],axis = 0)
            mesoscale_stats_df = mesoscale_stats_df.reindex(sorted(mesoscale_stats_df.columns), axis=1)

        self.targets_per_neuron   = targets_per_neuron
        self.somata               = somata

        return mesoscale_stats_df

    def asDict(self):
        return dict(
          somata = self.somata,
          targets_per_neuron = self.targets_per_neuron,
          neurite_length = self.neurite_length
        )

    def Done(self):
        self.somata = None
        self.targets_per_neuron = None
        self.neurite_length = None

    def barrel_specific_matrix(self, source_areas, target_areas, feature = 'length', extract = 'terminals'):


        barrel_stats = {}
        targets_per_neuron = OrderedDict()
        somata  = []

        acr_to_morpho_id = deepcopy(self.acr_to_morpho_id)

        annot_shape = self.annotation.shape
        out_orientation = ['um({})'.format(self.res),'PIR','corner']
        if self.neurite_length is not None:
            nld_list = self.neurite_length.keys()

        if 'L' in target_areas[0] or 'l' in target_areas[0] or 'layer' in target_areas[0]:
            barrel_annot, id2barrlabel = load_barrel_files(self.data_path, True)
        else:
            barrel_annot, id2barrlabel = load_barrel_files(self.data_path, False)

        for source_area in source_areas:

            if source_area not in acr_to_morpho_id:
                acr_to_morpho_id[source_area] = []
            for rand_area in acr_to_morpho_id.keys():  # First ensure that you have aggregated neurons of every child node
                if rand_area == '[background]' or rand_area == source_area: continue
                nu_id = str(self.acr2id[rand_area])
                if self.acr2id[source_area] in self.ancestorsById[nu_id]: # found a child area ...
                      acr_to_morpho_id[source_area].extend(acr_to_morpho_id[rand_area])
            if acr_to_morpho_id[source_area] == []: continue # nothing found here ...

            for neuron_id, neuron_soma_cord in acr_to_morpho_id[source_area]:
                # for each morphology, extract useful properties ...
                neuron_id = neuron_id.split('.')[0]
                nld_id = braintell_2_nld(nld_list, neuron_id)
                if nld_id == -1:
                    continue
                if self.un_num(self.neurite_length[nld_id]['soma']['region']) != self.un_num(self.neurite_length[nld_id]['soma']['correctedRegion']):
                    continue
                target_matches = [val for val in self.neurite_length[nld_id]['axon'].keys() if 'SSp' in val]
                if len(target_matches) == 0: # does not target any of the areas
                    continue

                db = 'mouselight' if 'AA' in neuron_id else 'braintell'
                neuronPath = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}.json.gz'.format(db,neuron_id)
                file_content = requests.get(neuronPath).content
                target_neuron = json.loads(zlib.decompress(file_content, 16+zlib.MAX_WBITS))

                neuron_cls = NeuronMorphology(neuronDict = target_neuron)
                neuron_name = neuron_id.split('.')[0]
                if db == 'braintell':
                    re_str = neuron_name.split('_')[1:3]
                    neuron_name = re_str[0]+'_'+ re_str[1]
                neuron_cls.transform(out_orientation)

                if feature == 'length':
                    if extract == 'terminals':
                        term = True
                    else:
                        term = False
                    axonal_length_barr = neuron_cls.compute_axonal_length(barrel_annot, id2barrlabel, terminals = True)
                    for key,val in axonal_length_barr.items():
                        if key in target_areas:
                            if (source_area,neuron_name) not in barrel_stats.keys():
                                barrel_stats[(source_area,neuron_name)] = OrderedDict()
                            if key not in barrel_stats[(source_area,neuron_name)].keys():
                                barrel_stats[(source_area,neuron_name)][key] = np.sum(val)
                            else:
                                barrel_stats[(source_area,neuron_name)][key] += np.sum(val)
                elif feature == 'counts':
                    anatomical_stats = neuron_cls.get_anatomical_stats(self.annotation,self.id2acr, mode = extract)
                    for key,val in anatomical_stats.items():
                        if key in target_areas:
                            if (source_area,neuron_name) not in barrel_stats.keys():
                                barrel_stats[(source_area,neuron_name)] = OrderedDict()
                            if key not in barrel_stats[(source_area,neuron_name)].keys():
                                barrel_stats[(source_area,neuron_name)][key] = len(val)
                            else:
                                barrel_stats[(source_area,neuron_name)][key] += len(val)

        barrel_stats_df = pd.DataFrame(barrel_stats).T
        if len(list(barrel_stats.keys())) > 0:
            barrel_stats_df.index = barrel_stats_df.index.rename(['Source','Neuron Id'])
            re_sort_cols = np.argsort([col[0] for col in barrel_stats_df.columns])
            re_sort_rows = np.argsort([col[0] for col in barrel_stats_df.index])
            barrel_stats_df = barrel_stats_df.reindex(barrel_stats_df.columns[re_sort_cols], axis = 1)
            barrel_stats_df = barrel_stats_df.reindex(barrel_stats_df.index[re_sort_rows],axis = 0)
            barrel_stats_df = barrel_stats_df.reindex(sorted(barrel_stats_df.columns), axis = 1)

        return barrel_stats_df
