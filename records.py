import itertools
import blastn


def get_hsp_features(record):
    "Return all of record's features with HSPs."
    return [f for f in record.features if hasattr(f, "hsp")]

def get_hsp_features_from(record, target):
    "Return all of record's features with HSPs that are aligned with target."
    return [f for f in get_hsp_features(record) if f.hsp.target.name == target.name]

def get_hsps(record):
    "Return all of record's HSPs."
    return [f.hsp for f in get_hsp_features(record)]

def get_hsps_from(record, target):
    "Return all of record's HSPs that are aligned with target."
    return [h for h in get_hsps(record) if h.target.name == target.name]

def has_hsp_from(record, target):
    "Return True if record has an HSP that is aligned with target, False otherwise."
    return bool(get_hsps_from(record, target))

def hsp_query(queries, hsp):
    return next(filter(lambda x: x.id == hsp.query.id, queries.values()), None)

def hsp_target(targets, hsp):
    return next(filter(lambda x: x.name == hsp.target.name, targets.values()), None)

def process_hsp(records, hsp):
    query = hsp_query(records, hsp)
    query_feature = blastn.hsp_to_feature(hsp, target=False)
    query.features.append(query_feature)

    target = hsp_target(records, hsp)
    target_feature = blastn.hsp_to_feature(hsp, target=True)
    target.features.append(target_feature)

    logger.info(f"Added HSP between {query.id} and {target.name}")
    return None

def has_transposon(record):
    return getattr(record, "transposon", False)

def get_transposon(record):
    return next(filter(lambda f: getattr(f, "transposon", False), record.features), None)

