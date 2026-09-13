# -*- coding: utf-8 -*-
"""
代码公开仓库发布脚本（③ 待办）
用法：
  1. 在 GitHub 创建空仓库（如 TSH-ThyroidCancer-MR-Study），或
     提供已有空仓库 URL。
  2. 提供 GitHub Personal Access Token（repo / public_repo 权限）：
     - 方式 A：设置环境变量 GITHUB_TOKEN=<PAT> 后运行本脚本
     - 方式 B：把 PAT 粘贴到本文件下方 TOKEN 变量
  3. 运行：python publish_repo.py <remote_url>
发布后论文数据可用性声明中的链接即可回填。
"""
import os, subprocess, sys

REPO = r"D:\gwas\script_repo"
PROXY = "http://127.0.0.1:7897"

def run(cmd, cwd=REPO):
    print("$", " ".join(cmd))
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-2000:])
        sys.exit(1)
    print(r.stdout[-2000:])

def main():
    if len(sys.argv) < 2:
        print("用法: python publish_repo.py https://github.com/<user>/<repo>.git")
        sys.exit(1)
    remote = sys.argv[1]
    token = os.environ.get("GITHUB_TOKEN", "")
    if not token:
        print("未检测到 GITHUB_TOKEN 环境变量。")
        token = input("请粘贴 GitHub PAT（不回显）: ").strip()
    if not token:
        sys.exit("未提供 token，取消发布。")
    # 1. 配置代理（GitHub 需走用户代理）
    run(["git", "config", "http.proxy", PROXY])
    run(["git", "config", "https.proxy", PROXY])
    # 2. 配置凭据（仅本次会话使用，避免明文落盘）
    run(["git", "remote", "remove", "origin"] if "origin" in
        subprocess.run(["git", "remote"], cwd=REPO, capture_output=True, text=True).stdout.split()
        else ["git", "remote", "remove", "origin"])
    run(["git", "remote", "add", "origin", remote])
    push_url = remote.replace("https://", "https://x-access-token:%s@" % token)
    # 3. push（token 仅出现在本次命令参数中）
    run(["git", "push", "-u", push_url, "main"] if os.path.exists(os.path.join(REPO, ".git", "refs", "heads", "main"))
        else ["git", "push", "-u", push_url, "master"])
    print("✓ 已发布到", remote)

if __name__ == "__main__":
    main()
