from Bio.SeqFeature import SimpleLocation

import copy

def overlapping(a, b):
    "Return True if location a overlaps with location b (ignoring strand)."
    return (a.start in b) or (a.end - 1 in b) or (b.start in a) or (b.end - 1 in a)

def features_in(location, record):
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
            new_location = SimpleLocation(start, end, strand=strand)
            new_feature.location = new_location
            new_features.append(new_feature)
    return new_features

def truncated_features_in(location, record):
    "Returns the features in record overlapping location and truncated."
    features = features_in(location, record)
    return truncated_features(features, location)

def translate_feature(feature, from_location, to_location):
    feature, = truncated_features([feature], from_location)
    start = to_location.start + feature.location.start - from_location.start
    end = to_location.end + feature.location.end - from_location.end
    strand = feature.location.strand
    
    location = SimpleLocation(min(start, end), max(start, end), strand=strand)
    feature = copy.deepcopy(feature)
    feature.location = location
    return feature

def translated_features_in(from_location, to_location, record):
    features = truncated_features_in(from_location, record)
    return [translate_feature(f, from_location, to_location) for f in features]
