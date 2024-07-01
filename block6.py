import matplotlib.pyplot as plt

# Definition of the TreeNode class, representing nodes in the tree
class TreeNode:
    def __init__(self, name=None, branch_length=None):
        self.name = name
        self.branch_length = branch_length
        self.children = []
        self.parent = None


# Function to split the Newick string into subtrees at the top level
def split_subtrees(s):
    parts, balance, current = [], 0, []
    for char in s:
        if char == '(':
            balance += 1
        elif char == ')':
            balance -= 1
        if char == ',' and balance == 0:
            parts.append(''.join(current).strip())
            current = []
        else:
            current.append(char)
    if current:
        parts.append(''.join(current).strip())
    return parts

# Function to extract branch name and length from a given string
def extract_branch_info(branch_part):
    if ':' in branch_part:
        name, length_str = branch_part.split(':', 1)
        try:
            return name.strip(), float(length_str)
        except ValueError:
            return name.strip(), None
    return branch_part.strip(), None

# Function to parse a Newick formatted string into a tree of TreeNode objects
def parse_newick(newick):
    def parse_subtree(subtree_str, parent=None):
        if not subtree_str:
            return None
        if subtree_str[0] == '(':
            close_index = subtree_str.rfind(')')
            branch_part = subtree_str[close_index+1:].strip()
            name, branch_length = extract_branch_info(branch_part)
            subtree_str = subtree_str[1:close_index].strip()
            node = TreeNode(name, branch_length)
            node.parent = parent
            for child_str in split_subtrees(subtree_str):
                child_node = parse_subtree(child_str, node)
                node.children.append(child_node)
            return node
        else:
            name, branch_length = extract_branch_info(subtree_str)
            return TreeNode(name, branch_length)

    return parse_subtree(newick.strip(';'))  # Strip the ending semicolon if present

# Function to count the number of labels (named nodes) in the tree
def count_labels(node):
    if node is None:
        return 0
    count = 1 if node.name else 0
    for child in node.children:
        count += count_labels(child)
    return count

# Function to draw the tree with square branches using matplotlib
def draw_square_branches_tree(node, label_count, base_font_size=20, x=0, y=0, x_offset=2, y_offset=2, axes=None, level=0):
    if node is None:
        return

    own_axes = False
    if axes is None:
        fig, axes = plt.subplots(figsize=(20, 20), dpi=100)
        own_axes = True

    font_size = max(base_font_size - (label_count // 4) * 2, 4)
    label = node.name if node.name else ''
    label_length = len(label)
    label_x = x + label_length * 0.05
    label_y = y
    axes.text(label_x, label_y, label, ha='right', va='center', fontsize=font_size, bbox=dict(facecolor='white', edgecolor='none', pad=1))

    y_step = y_offset / (2 ** (level / 1.15))

    for i, child in enumerate(node.children):
        child_x = x + x_offset
        child_y = y - y_step * (len(node.children) - 1) / 2 + i * y_step
        axes.plot([x, x + x_offset / 2], [y, y], 'k-')
        axes.plot([x + x_offset / 2, x + x_offset / 2], [y, child_y], 'k-')
        axes.plot([x + x_offset / 2, child_x], [child_y, child_y], 'k-')

        if child.branch_length is not None:
            label_pos_x = (x + x_offset / 2 + child_x) / 2
            label_pos_y = child_y - 0.1  # Adjust the y position of the branch length label
            axes.text(label_pos_x, label_pos_y, f'{child.branch_length}', fontsize=font_size, ha='center', va='top')

        draw_square_branches_tree(child, label_count, base_font_size, child_x, child_y, x_offset, y_offset, axes, level + 1)

    if own_axes:
        axes.set_xlim(left=0)
        axes.axis('off')
        plt.show()


# Example usage with a given Newick string
newick_str = final_tree
label_count = len(sequences)
tree = parse_newick(newick_str)
draw_square_branches_tree(tree, label_count)
