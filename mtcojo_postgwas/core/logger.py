"""
logger.py - Centralized Logging for mtcojo_postgwas

Features:
  - Dual output: coloured stdout + plain log file
  - Timestamped entries with log level
  - Colored step banners, PASS/WARN/FAIL indicators
  - Pipeline start banner with ASCII art
"""

import logging
import sys
import os
from datetime import datetime

# ── ANSI colour codes ─────────────────────────────────────────────────────────
R   = "\033[0m"          # reset
B   = "\033[1m"          # bold
DIM = "\033[2m"          # dim
GRN = "\033[38;5;82m"   # bright green
YEL = "\033[38;5;220m"  # amber yellow
RED = "\033[38;5;196m"  # bright red
CYN = "\033[38;5;51m"   # bright cyan
MGN = "\033[38;5;201m"  # magenta
BLU = "\033[38;5;39m"   # sky blue
ORG = "\033[38;5;214m"  # orange
WHT = "\033[38;5;255m"  # bright white


LOGO = f"""\
{MGN}{B}
  ┌─────────────────────────────────────────────────────────────────┐
  │                                                                 │
  │  {CYN}███╗   ███╗████████╗ ██████╗ ██████╗      ██╗ ██████╗ {MGN}         │
  │  {CYN}████╗ ████║╚══██╔══╝██╔════╝██╔═══██╗     ██║██╔═══██╗{MGN}         │
  │  {CYN}██╔████╔██║   ██║   ██║     ██║   ██║     ██║██║   ██║{MGN}         │
  │  {CYN}██║╚██╔╝██║   ██║   ██║     ██║   ██║██   ██║██║   ██║{MGN}         │
  │  {CYN}██║ ╚═╝ ██║   ██║   ╚██████╗╚██████╔╝╚█████╔╝╚██████╔╝{MGN}         │
  │  {CYN}╚═╝     ╚═╝   ╚═╝    ╚═════╝ ╚═════╝  ╚════╝  ╚═════╝ {MGN}         │
  │                                                                 │
  │  {WHT}mtcojo_postgwas{MGN} · GCTA mtCOJO Pipeline · PostGWAS Harmonisation │
  └─────────────────────────────────────────────────────────────────┘
{R}"""


# ── Formatter: coloured stdout ────────────────────────────────────────────────
class _ColorFormatter(logging.Formatter):
    COLORS = {
        logging.DEBUG:    DIM,
        logging.INFO:     WHT,
        logging.WARNING:  YEL,
        logging.ERROR:    RED,
        logging.CRITICAL: RED + B,
    }
    def format(self, record):
        c = self.COLORS.get(record.levelno, R)
        ts = datetime.now().strftime("%H:%M:%S")
        lvl = f"{record.levelname:<7}"
        return f"{DIM}{ts}{R} {c}{lvl}{R} {record.getMessage()}"


# ── Formatter: plain log file ─────────────────────────────────────────────────
class _FileFormatter(logging.Formatter):
    def format(self, record):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        return f"{ts} [{record.levelname:<8}] {record.getMessage()}"


def setup_logger(log_file: str = None, name: str = "mtcojo_postgwas") -> logging.Logger:
    """
    Create the pipeline logger.
    Attaches a coloured StreamHandler (stdout) and optionally a FileHandler.

    Args:
        log_file : Absolute path to the output .log file. Created automatically.
        name     : Logger name (default: 'mtcojo_postgwas').
    """
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)

    # Stdout handler (coloured)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(_ColorFormatter())
    logger.addHandler(sh)

    # File handler (plain text, always DEBUG level)
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(_FileFormatter())
        logger.addHandler(fh)
        logger.debug(f"Log file: {log_file}")

    logger.propagate = False
    return logger


def get_logger(name: str = "mtcojo_postgwas") -> logging.Logger:
    """Return existing logger (must call setup_logger first)."""
    return logging.getLogger(name)


# ── Convenience helpers ───────────────────────────────────────────────────────

def print_logo():
    print(LOGO)


def step_banner(logger: logging.Logger, title: str, step: int = None, total: int = None) -> None:
    step_str = f"[{step}/{total}] " if step and total else ""
    bar = "─" * 64
    logger.info(f"\n{CYN}{B}{bar}")
    logger.info(f"  {ORG}{B}{step_str}{WHT}{B}{title}")
    logger.info(f"{CYN}{B}{bar}{R}")


def log_pass(logger: logging.Logger, msg: str) -> None:
    logger.info(f"  {GRN}✔{R}  {msg}")


def log_warn(logger: logging.Logger, msg: str) -> None:
    logger.warning(f"  {YEL}⚠{R}  {msg}")


def log_fail(logger: logging.Logger, msg: str) -> None:
    logger.error(f"  {RED}✘{R}  {msg}")


def log_info(logger: logging.Logger, msg: str) -> None:
    logger.info(f"  {BLU}›{R}  {msg}")


def abort(logger: logging.Logger, msg: str) -> None:
    """Log a fatal error message and exit the pipeline."""
    logger.error("")
    for line in msg.strip().split("\n"):
        log_fail(logger, line)
    logger.error(f"\n  {RED}{B}━━━  PIPELINE ABORTED  ━━━{R}")
    raise SystemExit(1)


def log_cmd_script(out_prefix: str, step_title: str, cmd_str: str) -> None:
    """
    Appends an executed command to a single master shell script <out_prefix>_commands.sh.
    """
    cmd_file = f"{out_prefix}_commands.sh"
    os.makedirs(os.path.dirname(cmd_file), exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header_exists = os.path.exists(cmd_file)
    with open(cmd_file, "a", encoding="utf-8") as f:
        if not header_exists:
            f.write("#!/usr/bin/env bash\n")
            f.write("# " + "="*76 + "\n")
            f.write(f"# GCTA mtCOJO & LDSC Pipeline Executed Commands Log\n")
            f.write(f"# Generated: {ts}\n")
            f.write("# " + "="*76 + "\n\n")
        f.write(f"# ─── {step_title} [{ts}] ───\n")
        f.write(f"{cmd_str.strip()}\n\n")


def summary_table(logger: logging.Logger, rows: list, title: str = "Summary") -> None:
    """
    Print a boxed summary table.
    rows: list of (label, value) tuples.
    """
    label_w = max(len(r[0]) for r in rows) + 2
    val_w   = 36
    width   = label_w + val_w + 3
    logger.info(f"\n  {GRN}{B}╔{'═' * width}╗")
    logger.info(f"  {GRN}║  {WHT}{B}{title:<{width - 2}}{GRN}║")
    logger.info(f"  {GRN}╠{'═' * width}╣")
    for label, val in rows:
        label_str = f"  {label}".ljust(label_w)
        val_str   = str(val)
        if len(val_str) > val_w:
            val_str = "…" + val_str[-(val_w - 1):]
        logger.info(f"  {GRN}║{WHT}{label_str}{GRN}│ {WHT}{val_str:<{val_w}}{GRN}║")
    logger.info(f"  {GRN}╚{'═' * width}╝{R}")



