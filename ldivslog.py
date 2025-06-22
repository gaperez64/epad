import sys
import re
import yaml

start = re.compile(r"\[lindivs log\] start sys (\d+)")
end = re.compile(r"\[lindivs log\] end sys")
parent = re.compile(r"\[lindivs log\] from (\d+) (.*)")


class System:
    def __init__(self, idx):
        self.idx = idx
        self.reason = "given"
        self.system = []
        self.subsys = []


def read_log(fname):
    with open(fname, "r") as logf:
        root = None
        systems = dict()
        insys = False
        cur = None
        for line in logf:
            line = line.strip()
            match = start.match(line)
            if match:
                assert not insys
                insys = True
                idx = int(match.group(1))
                print(f"idx = {idx}")
                cur = System(idx)
                if root is None:
                    root = cur
                systems[idx] = cur
                continue

            match = parent.match(line)
            if match:
                assert insys
                pred = int(match.group(1))
                print(f"pred = {pred}")
                systems[pred].subsys.append(cur)
                cur.reason = match.group(2)
                continue

            match = end.match(line)
            if match:
                assert insys
                cur = None
                insys = False
                continue

            if line == "--":
                break

            # Default case
            cur.system.append(line)
    return root


if __name__ == "__main__":
    assert len(sys.argv) == 2
    tree = read_log(sys.argv[1])
    print(yaml.dump(tree, sort_keys=False))
