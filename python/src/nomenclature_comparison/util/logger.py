'''
Created on May 13, 2026

@author: pleyte
'''
import logging
import sys

def setup_logging(level=logging.INFO):
    """Configures the root logger for the entire nomenclature_comparison namespace."""
    logger = logging.getLogger("nomenclature_comparison")
    logger.setLevel(level)

    # Avoid adding multiple handlers if setup_logging is called twice
    if not logger.handlers:
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s - %(name)s - %(message)s'
        )
        
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger