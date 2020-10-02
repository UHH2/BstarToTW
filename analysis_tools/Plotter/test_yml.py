import yaml


config_file = "template_config.yml"
with open(config_file, 'r') as stream:
    try:
        print(yaml.safe_load(stream))
    except yaml.YAMLError as exc:
        print(exc)

