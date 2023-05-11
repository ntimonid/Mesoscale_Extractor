from cfg import *


def load_barrel_files(data_path, islayer = False):

    if islayer == True:
        barrel_annot, barr_hdr = nrrd.read(f'{data_path}/barrels_annotation_10_extended_layers.nrrd')
        id2barrlabel = pk.load(open(f'{data_path}/id2barrlabel_extended.pkl','rb'))
    else:
        barrel_annot, barr_hdr = nrrd.read(f'{data_path}/barrels_annotation_10_extended.nrrd')
        id2barrlabel = pk.load(open(f'{data_path}/id2barrlabel.json','rb'))

    return barrel_annot, id2barrlabel

def load_useful_variables(data_path, resolution = 10):

    try:

        [annotation,allenMeta] = nrrd.read('{}/annotation_{}.nrrd'.format(data_path,resolution))
        with open(os.path.join(data_path, 'ccf3_acr2id.json')) as fp:
             acr2id = json.load(fp)
        with open(os.path.join(data_path,'ancestorsById.json')) as fp:
            ancestorsById = json.load(fp)
    except:
        reference_space_key = 'annotation/ccf_2017'
        rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
        annotation, allenMeta = rspc.get_annotation_volume()
        tree = rspc.get_structure_tree(structure_graph_id=1)
        acr2id = tree.get_id_acronym_map()
        ancestorsById = tree.get_ancestor_id_map()

    id2acr = { id:acr for acr,id in acr2id.items() }
    if 0 not in id2acr:
      acr2id['[background]'] = 0
      id2acr[0] = '[background]'

    with open(os.path.join(data_path,'acr_to_morpho_id_new.pkl'), 'rb') as infile:
        acr_to_morpho_id = pk.load(infile)

    neurite_length = {}
    databases = ['braintell','mouselight']
    for dbName in databases:
          with open(os.path.join(data_path, 'neuriteLengthDistribution({}).json'.format(dbName))) as fp:
                neurite_length.update(json.load(fp))


    return annotation, neurite_length, acr_to_morpho_id, ancestorsById, id2acr


def braintell_2_nld(nld_list, morpho_id):
    if 'AA' in morpho_id:
        return morpho_id
    else:
        try:
            tmp = morpho_id.split('_reg')[0]
            matches =  [key for key in nld_list if tmp.split('_')[2] in key]
            return matches[0]
        except:
            return -1

def DownloadAxons(experimentId = None, mode = 'mouselight', rgb = '#FFBB00'): # Needs an update

    if mode == 'streamlines':
        save_path = os.path.join(main_path, 'Backup_Code/streamlines_tmp')
        infile = 'streamlines_{}.json.gz'.format(experimentId)
        streamline_link = 'https://neuroinformatics.nl/HBP/allen-connectivity-viewer/json/{}'.format(infile)
    elif mode == 'mouselight':
        save_path = os.path.join(main_path, 'Data Repositories/Mouselight/json')
        if type(experimentId) == int: # id not given directly
            experimentId += 1
            infile = 'AA{:04d}'.format(experimentId)
            experimentId = infile
        if '.json' not in experimentId:
            infile = '{}.json.gz'.format(experimentId)
        else:
            infile = experimentId
        streamline_link = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}'.format(mode,infile)
    elif mode == 'braintell':
        save_path = os.path.join(main_path, 'Data Repositories/Braintell')
        braintell_link = 'https://neuroinformatics.nl/HBP/braintell-viewer/fname2soma.bas(sba.ABA_v3(mm,RAS,ac)).json'
        response = urllib.request.urlopen(braintell_link)
        file_content = response.read()
        out = json.loads(file_content)
        if type(experimentId) == int: # id not given directly
            cnt = 0
            for key in out.keys():
                if cnt == experimentId:
                    experimentId = key
                    break
                cnt+=1
            infile = '{}.json.gz'.format(experimentId)
        if '.json.gz' not in experimentId:
            infile = '{}.json.gz'.format(experimentId)
        else:
            infile = experimentId
        streamline_link = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}'.format(mode,infile)

    infile2 = os.path.join(save_path,infile)
    print(infile, infile2, streamline_link)
    try:
        if os.path.exists(infile2) is False:
            wget.download(streamline_link,infile2)
        if '.json.gz' in infile2:
            with gzip.open(infile2, 'rb') as fp:
                file_content = fp.read()
        elif '.json' in infile2:
            with open(infile2, 'r') as fp:
                file_content = fp.read()
        out = json.loads(file_content)
    except:
        out = -1

    return out


def get_neuriteLengthDistribution():

    def lengthDistribution(neuron,neuriteType='axon',db='braintell',unit=1):
          typeFilter = { 'axon': (2,), 'dendrite': (3,4) }[neuriteType]
          points = neuron['treePoints']['data']
          lines = neuron['treeLines']['data']

          lengthDistr = {'@TOTAL': 0}
          oob = False
          for lineId,line in enumerate(lines):
            lineType,firstPoint,numPoints,prevLineId,negOffset = line
            if lineType in typeFilter:
              prevPoint_PIR = None
              if prevLineId:
                prevLine = lines[line[3]]
                prevPoint = points[prevLine[1]+prevLine[2]-1-line[4]]
                prevPoint_PIR = convert_PIR(prevPoint,db,unit)
              for p in range(line[1],line[1]+line[2]):
                point_PIR = convert_PIR(points[p],db,unit)
                point_10_PIR = np.round(point_PIR/10).astype(np.uint16)
                if np.any(point_10_PIR[0:3]>=ccf3_extent_10_PIR):
                  if not oob:
                    print('Out of bounds {}'.format(point_10_PIR))
                  oob = True # suppress further messages
                  continue
                id = annotationVolume[point_10_PIR[0],point_10_PIR[1],point_10_PIR[2]]
                if id not in id2acr:
                  print('Unknown id {}'.format(id))
                  continue

                if prevPoint_PIR is not None:
                  acr = id2acr[id]
                  if acr not in lengthDistr:
                    lengthDistr[acr] = 0
                  length_um = np.linalg.norm(point_PIR - prevPoint_PIR)
                  lengthDistr[acr] += length_um
                  lengthDistr['@TOTAL'] += length_um

                prevPoint_PIR = point_PIR

          return { k:int(v) for k,v in lengthDistr.items() }


    with gzip.open(fpath) as fp:
        neuron = json.loads(fp.read())
    # get soma coordinate (native)
    somaCoord = findSoma(neuron)
    unit = guessUnit(somaCoord)
    # transform to PIR, micrometers
    somaCoord_10_PIR = convert_PIR(somaCoord,db,unit,10)
    if np.any(somaCoord_10_PIR[0:3]>=ccf3_extent_10_PIR):
      print('Soma out of bounds {}'.format(somaCoord_10_PIR))
      somaRegion = '[?]'
    else:
      areaId = annotationVolume[somaCoord_10_PIR[0],somaCoord_10_PIR[1],somaCoord_10_PIR[2]]
      somaRegion = id2acr[areaId]

    region = id2acr[areaId]
    if db == 'braintell':
      # according to moesm4
      row = moesm4.loc[fid]
      correctedRegion = row['Manually_corrected_soma_region']
      correctedLayer = str(row['Cortical_layer'])
      if len(correctedLayer) and correctedLayer != 'nan':
        if correctedLayer == '6':
          correctedLayer = '6a'
          if region in regionsByLayer['Layer6b']: correctedLayer = '6b'
        if correctedLayer == '2/3':
          correctedLayer = '2_3'
        correctedRegion = areaPlusLayer(correctedRegion,correctedLayer)
    else:
      correctedRegion = region
      #correctedLayer = None
    neuriteLengthDistribution[fid] = {
      'fname': fname,
      'soma': {
        'region': somaRegion,
        'correctedRegion': correctedRegion,
        'corticalLayer': correctedLayer
      },
      'axon': lengthDistribution(neuron,'axon',db,unit),
      'dendrite': lengthDistribution(neuron,'dendrite',db,unit)
    }

    return neuriteLengthDistribution
