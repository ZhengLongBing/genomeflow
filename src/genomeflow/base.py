"""
序列的抽象基类。

设计理念：
- 提取 DNA 和 RNA 的共同行为到基类
- 使用 Python 的 ABC（Abstract Base Class）机制
- 子类只需实现特定的碱基规则和操作
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator


class InvalidSequenceError(ValueError):
    """当序列包含无效字符时抛出此异常。"""

    pass


class BaseSequence(ABC):
    """
    核酸序列的抽象基类。

    这是一个抽象类，不能直接实例化。
    子类必须实现 VALID_BASES 和 COMPLEMENT_MAP 属性。

    为什么使用抽象基类？
    1. 强制子类实现必要的方法
    2. 共享通用实现，减少代码重复
    3. 提供统一的接口，便于多态使用
    """

    @property
    @abstractmethod
    def VALID_BASES(self) -> frozenset[str]:
        """
        有效的碱基集合。

        子类必须实现此属性。
        DNA 返回 frozenset("ATGC")
        RNA 返回 frozenset("AUGC")
        """
        ...

    @property
    @abstractmethod
    def COMPLEMENT_MAP(self) -> dict[str, str]:
        """
        碱基配对映射。

        子类必须实现此属性。
        DNA: A↔T, G↔C
        RNA: A↔U, G↔C
        """
        ...

    def __init__(self, sequence: str) -> None:
        """
        创建序列对象。

        Args:
            sequence: 序列字符串（不区分大小写）

        Raises:
            InvalidSequenceError: 如果序列包含无效字符
        """
        normalized = sequence.upper()
        self._validate(normalized)
        self._sequence = normalized

    def _validate(self, sequence: str) -> None:
        """验证序列是否只包含有效碱基。"""
        invalid_chars = set(sequence) - self.VALID_BASES
        if invalid_chars:
            raise InvalidSequenceError(
                f"序列包含无效字符: {invalid_chars}。" f"有效字符: {self.VALID_BASES}"
            )

    @property
    def sequence(self) -> str:
        """返回原始序列字符串。"""
        return self._sequence

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index: int | slice) -> str:
        return self._sequence[index]

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequence)

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        if len(self._sequence) > 50:
            display = f"{self._sequence[:25]}...{self._sequence[-25:]}"
        else:
            display = self._sequence
        return f"{class_name}('{display}')"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, BaseSequence):
            return type(self) is type(other) and self._sequence == other._sequence
        if isinstance(other, str):
            return self._sequence == other.upper()
        return NotImplemented

    def __hash__(self) -> int:
        return hash((type(self).__name__, self._sequence))

    def complement(self) -> BaseSequence:
        """
        返回互补链。

        Returns:
            同类型的互补序列
        """
        comp_seq = "".join(self.COMPLEMENT_MAP[base] for base in self._sequence)
        return type(self)(comp_seq)

    def reverse_complement(self) -> BaseSequence:
        """
        返回反向互补链。

        Returns:
            同类型的反向互补序列
        """
        comp_seq = "".join(self.COMPLEMENT_MAP[base] for base in self._sequence)
        return type(self)(comp_seq[::-1])

    def gc_content(self) -> float:
        """
        计算 GC 含量。

        GC 含量是 G 和 C 碱基的百分比，对 DNA 和 RNA 都有意义。

        Returns:
            GC 含量，范围 0.0-1.0
        """
        if len(self._sequence) == 0:
            return 0.0
        gc_count = self._sequence.count("G") + self._sequence.count("C")
        return gc_count / len(self._sequence)
