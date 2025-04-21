import copy
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

def runners(a):
    sublists = []
    i = 0
    while i < len(a):
        sublist = [a[i]]
        while i < len(a) - 1 and a[i] + 1 == a[i + 1]:
            i = 1 + i
            sublist.append(a[i])
        sublists.append(sublist)
        i = 1 + i

    return sublists

def unannotated_locations(record):
    missing = []
    for i in range(len(record)):
        if not any(i in feature for feature in record.features):
            missing.append(i)

    sections = runners(missing)
    locations = []
    for section in sections:
        location = FeatureLocation(section[0], section[-1] + 1)
        locations.append(location)
    return locations

def overlapping(a, b):
    "Return True if location a overlaps with location b (ignoring strand)."
    return (a.start in b) or (a.end - 1 in b) or (b.start in a) or (b.end - 1 in a)

def features_in(record, location):
    "Returns the features in record that overlap with location."
    return [f for f in record.features if overlapping(location, f.location)]

def truncated_features(features, location):
    "Truncates the locations of features so that they are within location."
    new_features = []
    for feature in features:
        if overlapping(feature.location, location):
            start = max(feature.location.start, location.start)
            end = min(feature.location.end, location.end)
            strand = feature.location.strand
            new_feature = copy.deepcopy(feature)
            new_location = FeatureLocation(start, end, strand=strand)
            new_feature.location = new_location
            new_features.append(new_feature)
    return new_features
    
def region_between(record, x, y):
    "Returns a new record substituting N at positions i for i < x or i >= y."
    record = copy.deepcopy(record)
    old_seq = list(record.seq)
    new_seq = ["N" if x > i or i >= y else old_seq[i] for i in range(len(old_seq))]
    record.seq = Seq("".join(new_seq))
    return record

def translate_feature(feature, from_location, to_location):
    feature, = truncated_features([feature], from_location)
    start = to_location.start + feature.location.start - from_location.start
    end = to_location.end + feature.location.end - from_location.end
    strand = feature.location.strand
    
    location = SimpleLocation(min(start, end), max(start, end), strand=strand)
    feature = copy.deepcopy(feature)
    feature.location = location
    return feature
