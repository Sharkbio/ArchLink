import logging
import sys
import os
from contrastive_learning import contrastive_learning_main
from generate01 import generate_main
from binning02 import binning_main
from generate_cos03 import generate_cos_main
from connect04 import connect_main
import argparse
import yaml
import multiprocessing
import shutil
from pathlib import Path


def load_yaml_with_vars(path):
    """load YAML file and process variable substitutions like {var}"""
    with open(path, "r") as f:
        data = yaml.safe_load(f)
    global_env = {}
    def collect_env(d):
        if isinstance(d, dict):
            for k, v in d.items():
                if isinstance(v, (str, int, float)):
                    global_env[k] = str(v)
                collect_env(v)
        elif isinstance(d, list):
            for item in d:
                collect_env(item)

    collect_env(data)

    def expand(value):
        if isinstance(value, str):
            for k, v in global_env.items():
                value = value.replace(f"{{{k}}}", v)
            return value
        return value

    def recursive_expand(d):
        if isinstance(d, dict):
            return {k: recursive_expand(expand(v)) for k, v in d.items()}
        elif isinstance(d, list):
            return [recursive_expand(expand(i)) for i in d]
        else:
            return expand(d)

    return recursive_expand(data)

class Args:
    """
    auto initialize from YAML
    """

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

        if hasattr(self, "output_path"):
            self.data = self.output_path

    @classmethod
    def from_yaml(cls, yaml_path):
        cfg = load_yaml_with_vars(yaml_path)

        #  common/path
        path_cfg   = cfg.get("common", {}).get("path", {})
        share_cfg  = cfg.get("contrastive_learning", {}).get("share_params", {})
        train_cfg  = cfg.get("contrastive_learning", {}).get("train", {})
        linking_binning_cfg  = cfg.get("linking", {}).get("binning", {})
        trans_cfg  = cfg.get("linking", {}).get("prodigal", {})
        model_cfg  = cfg.get("linking", {}).get("model", {})
        
        


        merged = {}
        merged.update(path_cfg)
        merged.update(share_cfg)
        merged.update(train_cfg)
        merged.update(linking_binning_cfg)
        merged.update(trans_cfg)
        merged.update(model_cfg)

        return cls(**merged)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run the ArchLink metagenomic binning and linking pipeline."
    )
    parser.add_argument(
        "--config",
        default="configuration.yaml",
        help="Path to the YAML configuration file. Defaults to configuration.yaml in the repository root.",
    )
    return parser


def setup_logging(logger_name, output_path):
    """set logging to file and console"""
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(message)s')
    
    # 控制台输出
    console_hdr = logging.StreamHandler(sys.stdout)
    console_hdr.setFormatter(formatter)
    # logger.addHandler(console_hdr)

    # 文件输出
    try:
        os.makedirs(output_path, exist_ok=True)
        log_filepath = os.path.join(output_path, 'ArchLink.log')
        handler = logging.FileHandler(log_filepath)
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.info(f"File logging enabled at: {log_filepath}")
    except Exception as e:
        logger.warning(f"Could not set up file logging: {e}")
        
    return logger


def remove_temp_file(file_path):
    # 删除 aug 文件夹
    aug_path = os.path.join(file_path, 'aug')
    if os.path.exists(aug_path):
        shutil.rmtree(aug_path)  

    # 删除 cluster_res 文件夹
    cluster_res_path = os.path.join(file_path, 'cluster_res')
    if os.path.exists(cluster_res_path):
        shutil.rmtree(cluster_res_path)

    # 删除 contrastive_learning_res_bins 文件夹
    contrastive_learning_res_path = os.path.join(file_path, 'contrastive_learning_res_bins')
    if os.path.exists(contrastive_learning_res_path):
        shutil.rmtree(contrastive_learning_res_path)


def main():
    try:
        multiprocessing.set_start_method("spawn", force=True)
    except RuntimeError as e:
        # 忽略已经设置的错误
        if "context has already been set" not in str(e):
            print(f"warning：can't set multiprocessing start method: {e}")
        pass
    
    cli_args = build_parser().parse_args()
    config_path = Path(cli_args.config).expanduser().resolve()
    args = Args.from_yaml(str(config_path))
    
    # 初始化日志
    logger = setup_logging('ArchLink', args.output_path)
    logger.info("ArchLink starting full pipeline.")
    logger.info(f"Configuration file: {config_path}")
    logger.info(f"Output Directory: {args.output_path}")
    contrastive_learning_main.contrastive_learning_main(args,logger)
    logger.info("STEP:generate...")
    generate_main.generate_init(args)
    logger.info("binning ...")
    binning_main.binning_init(args,logger)
    logger.info("generate cos ...")
    generate_cos_main.main(args)
    logger.info("connect ...")
    connect_main.connect_main(args)
    
    # remove_temp_file(args.output_path)


if __name__ == '__main__':
    main()
