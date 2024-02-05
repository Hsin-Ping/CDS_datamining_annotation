import logging
import os

logging.basicConfig(level=logging.INFO, filename="log.out", format="%(asctime)s - %(name)s %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
