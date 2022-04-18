import yaml
NCBI_CFG_FILE = "../CONFIG_NCBI.yml"

NCBI_CFGS = yaml.load(open(NCBI_CFG_FILE, "r"), Loader=yaml.BaseLoader)
